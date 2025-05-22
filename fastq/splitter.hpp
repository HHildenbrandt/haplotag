/* fastq/splitter.hpp
 *
 * Copyright (c) 2025 Hanno Hildenbrandt <h.hildenbrandt@rug.nl>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
*/

#pragma once

#include <cassert>
#include <string>
#include <array>
#include <string_view>
#include <type_traits>
#include <functional>
#include <memory>
#include <thread>
#include <bit>
#include <vector>
#include "reader.hpp"


namespace fastq {

  using recycle_fun = std::function<void(std::queue<chunk_t>&&)>;


  // return type of Splitter::operator()(size_t n).
  // shares ownership of the viewed memory (a.k.a. reader-chunks).
  template <typename T>
  class blk_reads_t {
  public:
    blk_reads_t(blk_reads_t&&) = default;
    blk_reads_t& operator=(blk_reads_t&&) = default;
    
    blk_reads_t(std::unique_ptr<T[]>&& val, size_t size, auto&& storage) : 
      views_(std::move(val)), 
      size_(size),
      shared_storage_(std::move(storage)) {
    }

    bool empty() const noexcept { return 0 == size_; }
    size_t size() const noexcept { return size_; }

    T operator[](size_t i) const { return views_[i]; }
    const T* data() const noexcept { return views_.get(); }
    const T* begin() const noexcept { return views_.get(); }
    const T* end() const noexcept { return views_.get() + size_; }

  private:
    std::unique_ptr<T[]> views_;
    const size_t size_ = 0;
    std::vector<chunk_t> shared_storage_;   // keep chunks alive
  };
  

  template <
    typename Reader,
    typename TrimPolicy,  // trim tail of new chunk
    typename SplitPolicy
  >
  class base_splitter {
  public:
    using value_type = SplitPolicy::value_type;
    using blk_type = blk_reads_t<value_type>;
    
    explicit base_splitter(const std::filesystem::path& path) : reader_(std::make_shared<Reader>(path)) {}

    bool eof() const noexcept { return last_ && cv_.empty(); }
    bool failed() const noexcept { return reader_->failed(); }
    const Reader& reader() const noexcept { return *reader_.get(); }

    // returns view into memory we don't own
    // valid until next call to operator()
    value_type operator()() {
      if (cv_.empty()) [[unlikely]] {
        if (!next_chunk()) return {};
      }
      return SplitPolicy::apply(cv_);
    }

    // returns views into up to n items
    // valid over the live time of the returned object
    blk_reads_t<value_type> operator()(size_t n) {
      using allocator = Reader::allocator;
      auto _ = buffer_guard();
      auto v = std::make_unique<value_type[]>(n);
      size_t i = 0;
      for (; !eof() && (i < n); ++i) {
        v[i] = this->operator()();
      }
      return blk_reads_t<value_type>{ std::move(v), i, release_buffer()};
    }

  private:
    bool next_chunk() {
      if (last_) [[unlikely]] {
        return false;
      }
      auto chunk = (*reader_.get())();
      last_ = chunk.last;
      if (!chunks_.empty()) {
        // copy over into window spare
        auto& prev = chunks_.back();
        std::memcpy(chunk.data() - tail_len_, prev.data() + prev.size - tail_len_, tail_len_);
        if (!buffered_) {
          chunks_.pop_back();
        }
      }
      auto tcv = TrimPolicy::apply(chunk);
      cv_ = { chunk.cv().data() - tail_len_, tcv.length() + tail_len_ };
      tail_len_ = chunk.cv().length() - tcv.length();
      assert(tail_len_ < Reader::window);
      chunks_.push_back(std::move(chunk));
      return true;
    }

    std::vector<chunk_t> release_buffer() {
      auto tmp = std::move(chunks_); chunks_ = {};
      chunks_.push_back(tmp.back());   // keep current
      return tmp;
    }

    struct buffer_guard_t{ bool& buffered; ~buffer_guard_t() { buffered = false; } };
    auto buffer_guard() { return buffer_guard_t(buffered_ = true); }

    std::string_view cv_;          // view into current chunk
    std::vector<chunk_t> chunks_;   // buffered chunks
    size_t tail_len_ = 0;
    bool last_ = false;
    bool buffered_ = false;
    std::shared_ptr<Reader> reader_;
  };


  namespace policy {

    struct char_trim {
      static std::string_view apply(const chunk_t& chunk) noexcept {
        return chunk.cv();
      }
    };


    struct char_apply {
      using value_type = char;

      static value_type apply(std::string_view& /* in/out */ cv) noexcept {
        auto ret = cv[0];
        cv.remove_prefix(1);
        return ret;
      }
    };


    template <
      char StartSymbol,
      char StopSymbol
    >
    struct brace_trim {
      static constexpr char delim[] = { StopSymbol, StartSymbol, '\0' };

      static std::string_view apply(const chunk_t& chunk) noexcept {
        auto cv = chunk.cv();
        if (!chunk.last) [[likely]] {
          // not the last chunk: 
          // identify tail section potentially overlapping into next chunk
          const auto p1 = cv.rfind(delim);
          if (p1 != cv.npos) [[likely]] {
            cv = cv.substr(0, p1 + 1);
          }
        }
        return cv;
      };
    };


    template <
      char StartSymbol,
      char StopSymbol,
      bool SkipStopSymbol
    >
    struct brace_split {
      static constexpr char delim[] = { StopSymbol, StartSymbol, '\0' };
      using value_type = std::string_view;

      static value_type apply(std::string_view& /* in/out */ cv) noexcept {
        assert(cv[0] == StartSymbol);
        auto ret = std::string_view{};
        const auto p1 = cv.find(delim);
        if (p1 == cv.npos) [[unlikely]] {
          ret = cv;
          cv = {};
        }
        else {
          ret = cv.substr(0, p1 + !SkipStopSymbol);
          cv.remove_prefix(p1 + 1);
        }
        return ret;
      }
    };


    template <
      char Delim,         // delimiter
      bool IsStartSymbol  // delimiter is start symbol
    >
    struct delim_trim {
      // sets cv and tail for a new chunk.
      static std::string_view apply(const chunk_t& chunk) noexcept {
        auto cv = chunk.cv();
        if (!chunk.last) [[likely]] {
          auto p1 = cv.rfind(Delim);
          cv = cv.substr(0, p1 + !IsStartSymbol);
        }
        return cv;
      };
    };


    template <
      char Delim,
      bool IsStartSymbol,
      bool SkipDelim
    >
    struct delim_split {
      using value_type = std::string_view;

      static value_type apply(std::string_view& /* in/out */ cv) noexcept {
        auto ret = value_type{};
        if constexpr (IsStartSymbol) {
          assert(cv[0] == Delim);
          const auto p1 = cv.substr(1).find(Delim);
          ret = cv.substr(SkipDelim, p1);
          cv = cv.substr(ret.length() + (ret.length() != cv.length()));
        }
        else {
          const auto p1 = cv.find(Delim);
          assert(p1 != cv.npos);
          ret = cv.substr(0, p1 + (1 - SkipDelim));
          cv = cv.substr(ret.length() + 1);
        }
        return ret;
      }
    };


    template <size_t Mask>
    struct masked_seq_split {
      static constexpr size_t mask = Mask;
      static constexpr int size = std::popcount(Mask);
      using value_type = std::array<std::string_view, size>;

      static value_type apply(std::string_view& /* in out */ cv) noexcept {
        using field_split = delim_split<'\n', false, true>;
        auto ret = value_type{};
        unsigned j = 0;
        for (auto i = 0; i < 4; ++i) {
          auto field = field_split::apply(cv);
          if (Mask & (1u << i)) {
            ret[j++] = field;
          }
        }
        return ret;
      }
    };

  }


  // the only copying splitter in our arsenal
  template <typename Reader>
  using char_splitter = base_splitter<
    Reader,
    policy::char_trim,
    policy::char_apply
  >;


  template <
    typename Reader,
    char StartSymbol,
    char StopSymbol,
    bool SkipStopSymbol
  >
  using brace_splitter = base_splitter<
    Reader,
    policy::brace_trim<StartSymbol, StopSymbol>,
    policy::brace_split<StartSymbol, StopSymbol, SkipStopSymbol>
  >;


  template <
    typename Reader,
    char Delim,
    bool IsStartSymbol,
    bool SkipDelim
  >
  using delim_splitter = base_splitter<
    Reader, 
    policy::delim_trim<Delim, IsStartSymbol>,
    policy::delim_split<Delim, IsStartSymbol, SkipDelim>
  >;
  
  
  template <typename Reader = reader_t>
  using line_splitter = delim_splitter<Reader, '\n', false, true>;
  

  template <typename Reader = reader_t>
  using seq_splitter = brace_splitter<Reader, '@', '\n', false>;
  

  // masked fasta sequence splitter
  // a '1' in the bitset 'Mask' selects the field to keep
  // e.g. Mask = 0b1111 keep all 4 fields
  //      Mask = 0b1110 drop first field, keep 2ndm 3rd and 4th
  template <size_t Mask = 0b1111, typename Reader = reader_t>
  using seq_field_splitter = base_splitter<
    Reader, 
    policy::brace_trim<'@', '\n'>,
    policy::masked_seq_split<Mask>
  >;

}
