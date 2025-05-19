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
#include <queue>
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
    std::queue<chunk_t> shared_storage_;   // keep chunks alive
  };
  

  template <
    typename Reader,
    typename TrimPolicy,
    typename SplitPolicy
  >
  class base_splitter {
  public:
    using value_type = SplitPolicy::value_type;
    using blk_type = blk_reads_t<value_type>;
    
    explicit base_splitter(const std::filesystem::path& path) : reader_(std::make_shared<Reader>(path)) {
      chunks.emplace();
    }

    bool eof() const noexcept { return chunks.back().last && cv.empty(); }
    bool failed() const noexcept { return reader_->failed(); }
    const Reader& reader() const noexcept { return *reader_.get(); }

    // returns view into memory we don't own
    // valid until next call to operator()
    value_type operator()() {
      if (cv.empty()) [[unlikely]] {
        if (!next_chunk()) return {};
      }
      return SplitPolicy::apply(cv);
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
      if (chunks.back().last) [[unlikely]] {
        return false;
      }
      auto& chunk = chunks.emplace((*reader_.get())(tail));
      TrimPolicy::apply(chunk, cv, tail);
      if (!buffered) {
        chunks.pop();
      }
      return true;
    }

    std::queue<chunk_t> release_buffer() {
      auto tmp = std::move(chunks); chunks = {};
      chunks.push(tmp.back());   // keep current
      return tmp;
    }

    struct buffer_guard_t{ bool& buffered; ~buffer_guard_t() { buffered = false; } };
    auto buffer_guard() { return buffer_guard_t(buffered = true); }

    std::string_view cv;       // view into current chunk
    std::string_view tail;     // orphan last chunk
    std::queue<chunk_t> chunks;   // buffered chunks
    bool buffered = false;

  private:
    std::shared_ptr<Reader> reader_;
  };


  namespace policy {

    struct char_trim {
      static void apply(const chunk_t& chunk, std::string_view& /* out */ cv, std::string_view& /* out */ tail) noexcept {
        cv = std::string_view(chunk);
        tail = {};
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
      char Delim,
      bool SkipDelim,
      bool StartSymbol
    >
    struct delim_trim {
      // sets cv and tail for a new chunk.
      static void apply(const chunk_t& chunk, std::string_view& /* out */ cv, std::string_view& /* out */ tail) noexcept {
        cv = std::string_view{chunk};
        tail = {};
        if (!chunk.last) [[likely]] {
          // not the last chunk: 
          // identify tail section potentially overlapping into next chunk
          if constexpr (StartSymbol) {
            auto p1 = cv.rfind(Delim);
            tail = cv.substr(p1);
            cv = cv.substr(0, p1);
          }
          else {
            // stop symbol
            auto p1 = cv.rfind(Delim);
            tail = cv.substr(p1 + 1);
            cv = cv.substr(0, p1 + 1);
          }
        }
      };
    };


    template <
      char Delim,
      bool SkipDelim,
      bool StartSymbol
    >
    struct delim_split {
      using value_type = std::string_view;

      static value_type apply(std::string_view& /* in/out */ cv) noexcept {
        auto ret = value_type{};
        if constexpr (StartSymbol) {
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

      static value_type apply(std::string_view& cv) noexcept {
        auto sv = delim_split<'@', false, true>::apply(cv);
        auto ret = value_type{};
        unsigned j = 0;
        for (auto i = 0; i < size; ++i) {
          auto p1 = sv.find('\n');
          if (Mask & (1u << i)) {
            ret[j++] = sv.substr(0, p1);
          }
          sv.remove_prefix(p1 + 1);
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
    char Delim,
    bool SkipDelim,
    bool StartSymbol,
    typename Reader = reader_t
  >
  using delim_splitter = base_splitter<
    Reader, 
    policy::delim_trim<Delim, SkipDelim, StartSymbol>,
    policy::delim_split<Delim, SkipDelim, StartSymbol>
  >;
  
  
  template <typename Reader = reader_t>
  using line_splitter = delim_splitter<'\n', true, false, Reader>;
  

  template <typename Reader = reader_t>
  using seq_splitter = delim_splitter<'@', false, true, Reader>;
  

  // masked fasta sequence splitter
  // a '1' in the bitset 'Mask' selects the field to keep
  // e.g. Mask = 0b1111 keep all 4 fields
  //      Mask = 0b1110 drop first field, keep 2ndm 3rd and 4th
  template <size_t Mask = 0b1111, typename Reader = reader_t>
  using seq_field_splitter = base_splitter<
    Reader, 
    policy::delim_trim<'@', false, true>,
    policy::masked_seq_split<Mask>
  >;
  
}
