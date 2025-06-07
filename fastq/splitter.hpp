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
#include <utility>
#include <string>
#include <array>
#include <string_view>
#include <type_traits>
#include <functional>
#include <memory>
#include <thread>
#include <bit>
#include <vector>
#include "fastq.hpp"
#include "reader.hpp"


namespace fastq {

  // return type of Splitter::operator()(size_t n).
  // shares ownership of the viewed memory (a.k.a. reader-chunks).
  template <typename T>
  class blk_reads_t {
  public:
    blk_reads_t(blk_reads_t&&) = default;
    blk_reads_t& operator=(blk_reads_t&&) = default;
    
    blk_reads_t(std::vector<T>&& val, auto&& storage) : 
      val_{std::move(val)}, 
      shared_storage_(std::move(storage)) {
    }

    bool empty() const noexcept { return val_.empty(); }
    size_t size() const noexcept { return val_.size(); }

    T operator[](size_t i) const { return val_[i]; }
    const T* data() const noexcept { return val_.data(); }
    auto begin() const noexcept { return val_.begin(); }
    auto end() const noexcept { return val_.end(); }

  private:
    std::vector<T> val_;
    std::vector<chunk_t> shared_storage_;   // keep chunks alive
  };
  

  template <
    typename TrimPolicy,  // trim tail of new chunk
    typename SplitPolicy
  >
  class chunk_splitter {
  public:
    using trim_policy = TrimPolicy;
    using split_policy = SplitPolicy;
    using value_type = split_policy::value_type;

    // assigns new chunk and returns old chunk
    chunk_t& assign(chunk_t&& chunk) {
      if (tail_len_) {
        // copy tail over into window spare
        std::memcpy(chunk.data() - tail_len_, chunk_.data() + chunk_.size - tail_len_, tail_len_);
      }
      auto trimmed_cv = trim_policy::apply(chunk);
      cv_ = { chunk.cv().data() - tail_len_, trimmed_cv.length() + tail_len_ };
      tail_len_ = chunk.cv().length() - trimmed_cv.length();
      chunk_ = std::move(chunk);
      return chunk_;
    }

    str_view cv() const noexcept { return cv_; }
    chunk_t chunk() const noexcept { return chunk_; }
    bool empty() const noexcept { return cv_.empty(); }

    value_type operator()() {
      assert(!empty());
      return split_policy::apply(cv_);
    }

  private:
    str_view cv_;   // view into chunk_
    chunk_t chunk_;         // current chunk
    size_t tail_len_ = 0;
  };


  template <
    typename Reader,
    typename ChunkSplitter
  >
  class base_splitter {
  public:
    using value_type = ChunkSplitter::value_type;
    using blk_type = blk_reads_t<value_type>;

    base_splitter() = default;
    base_splitter(base_splitter&&) = default;
    base_splitter& operator=(base_splitter&&) = default;

    explicit base_splitter(const std::filesystem::path& path) : reader_(std::make_shared<Reader>(path)) {}

    bool eof() const noexcept { return last_ && chunk_splitter_.empty(); }
    bool failed() const noexcept { return !reader_ || reader_->failed(); }
    const Reader& reader() const noexcept { return *reader_.get(); }

    // returns view into memory we don't own
    // valid until next call to operator()
    value_type operator()() {
      if (chunk_splitter_.empty()) [[unlikely]] {
        if (!next_chunk()) return {};
      }
      return chunk_splitter_();
    }

    // returns views into up to n items
    // valid over the live time of the returned object
    blk_reads_t<value_type> operator()(size_t n) {
      using allocator = Reader::allocator_t;
      auto v = std::vector<value_type>{};
      v.reserve(n);
      shared_storage_ = std::vector<chunk_t>{chunk_splitter_.chunk()};
      buffer_guard _{buffered_ = true};
      size_t i = 0;
      for (; !eof() && (i < n); ++i) {
        v.emplace_back(this->operator()());
      }
      return blk_reads_t<value_type>{ std::move(v), std::move(shared_storage_)};
    }

  private:
    bool next_chunk() {
      if (last_) [[unlikely]] {
        return false;
      }
      last_ = chunk_splitter_.assign((*reader_.get())()).last;
      if (buffered_) {
        shared_storage_.push_back(chunk_splitter_.chunk());
      }
      return true;
    }

    struct buffer_guard{ bool& buffered; ~buffer_guard() { buffered = false; } };

    ChunkSplitter chunk_splitter_;
    bool last_ = false;
    bool buffered_ = false;
    std::vector<chunk_t> shared_storage_;   // temporary keep alive buffer
    std::shared_ptr<Reader> reader_;
  };


  namespace policy {

    struct char_chunk_trim {
      static str_view apply(const chunk_t& chunk) noexcept {
        return chunk.cv();
      }
    };


    struct char_apply {
      using value_type = char;

      static value_type apply(str_view& /* in/out */ cv) noexcept {
        auto ret = cv[0];
        cv.remove_prefix(1);
        return ret;
      }
    };


    template <typename Chr, Chr Delim>
    struct delim_chunk_trim {
      static str_view apply(const chunk_t& chunk) noexcept {
        auto cv = chunk.cv();
        if (!chunk.last) [[likely]] {
          if (auto p1 = cv.rfind(Delim); p1 != cv.npos) [[likely]] {
            cv = cv.substr(0, p1 + 1);
          }
        }
        return cv;
      };
    };


    template <
      typename Chr, Chr Delim,
      int RemoveFront,      // # leading characters to remove
      int RemoveBack        // # trailing characters to remove
    >
    struct delim_split {
      using value_type = str_view;

      static value_type apply(str_view& /* in/out */ cv) noexcept {
        auto ret = str_view{};
        if (auto p1 = cv.find(Delim); p1 != cv.npos) [[likely]] {
          assert(cv.length() > (RemoveFront + RemoveBack));
          ret = cv.substr(RemoveFront, p1 + (1 - (RemoveFront + RemoveBack)));
          cv.remove_prefix(p1 + 1);
        } 
        else {
          std::swap(ret, cv);
        }
        return ret;
      }
    };


    template <size_t N, size_t Mask>
    struct masked_lines_split {
      static_assert(N < 64);
      static constexpr size_t mask = Mask;
      static constexpr int size = std::popcount(Mask);
      using value_type = std::array<str_view, size>;

      static value_type apply(str_view& /* in out */ cv) noexcept {
        using field_split = delim_split<char, '\n', 0, 1>;
        auto ret = value_type{};
        auto it = ret.begin();
        for (auto i = 0; i < N; ++i) {
          auto field = field_split::apply(cv);
          if (Mask & (1u << i)) {
            *it++ = field;
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
    chunk_splitter<
      policy::char_chunk_trim,
      policy::char_apply
    >
  >;


  template <typename Reader = reader_t>
  using line_splitter = base_splitter<
    Reader,
    chunk_splitter<
      policy::delim_chunk_trim<char, '\n'>,
      policy::delim_split<char, '\n', 0, 1>
    >
  >;
  
  // char* as template parameter better have
  // static linkage...
  static constexpr char SeqDelim[] = "\n@";

  template <typename Reader = reader_t>
  using seq_splitter = base_splitter<
    Reader,
    chunk_splitter<
      policy::delim_chunk_trim<const char*, SeqDelim>,
      policy::delim_split<const char*, SeqDelim, 0, 1>
    >
  >;
  

  // masked fasta sequence splitter
  // a '1' in the bitset 'Mask' selects the field to keep
  // e.g. Mask = 0b1111 keep all 4 fields
  //      Mask = 0b1110 drop first field, keep 2ndm 3rd and 4th
  template <size_t Mask = 0b1111, typename Reader = reader_t>
  using seq_field_splitter = base_splitter<
    Reader, 
    chunk_splitter<
      policy::delim_chunk_trim<const char*, SeqDelim>,
      policy::masked_lines_split<4, Mask>
    >
  >;


}
