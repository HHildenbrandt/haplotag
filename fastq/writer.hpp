/* fastq/writer.hpp
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

/*
 * Loosely based on Mark Adler's pigz code: 
 * https://github.com/madler/pigz.git
 * 
 * Pretty much everything zlib-related is rooted on Mark's work.
 * Thus, all credits to him!
*/

#pragma once

#include <cassert>
#include <cstring>
#include <stdexcept>
#include <filesystem>
#include <fstream>
#include <vector>
#include <span>
#include <string_view>
#include <thread>
#include <atomic>
#include <memory>
#include <utility>
#include <bit>
#include <device/pool.hpp>
#include <zlib.h>
#include "fastq.hpp"


namespace fastq {


  namespace {

    // initialize a raw deflate stream (no header, no checksum)
    inline void z_stream_init(z_stream& strm) {
      std::memset(&strm, 0, sizeof(z_stream));
      auto ret = deflateInit2(&strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED, -15, 9, Z_DEFAULT_STRATEGY);
      if (ret == Z_MEM_ERROR) {
        throw std::runtime_error("fastq::z_stream_init: not enough memory");
      }
      if (ret != Z_OK) {
        throw std::runtime_error("fastq::z_stream_init: internal error");
      }
    }


    inline void z_stream_reset(z_stream& strm) {
      (void)deflateReset(&strm);
      (void)deflateParams(&strm, Z_DEFAULT_COMPRESSION, Z_DEFAULT_STRATEGY);
    }


    // returns val in little endian byte order
    // inverse of the common 'hton' function
    // gz footer is little endian for some reason 
    constexpr uint32_t htogz(uint32_t val) {
      if constexpr (std::endian::native == std::endian::big) {
#if __cpp_lib_byteswap       
        return std::byteswap(val);
#else
        return (val >> 24) |
               ((val << 8) & 0x00FF0000) |
               ((val >> 8) & 0x0000FF00) |
               (val << 24);
#endif
      }
      return val;
    }


    // returns val in little endian byte order
    // inverse of the common 'hton' function
    // gz footer is little endian for some reason 
    constexpr uint64_t htogz(uint64_t val) {
      if constexpr (std::endian::native == std::endian::big) {
#if __cpp_lib_byteswap       
        return std::byteswap(val);
#else
        return (val >> 56) |
               ((val << 40) & 0x00FF000000000000) |
               ((val << 24) & 0x0000FF0000000000) |
               ((val << 8) & 0x000000FF00000000) |
               ((val >> 8) & 0x00000000FF000000) |
               ((val >> 24) & 0x0000000000FF0000) |
               ((val >> 40) & 0x000000000000FF00) |
               (val << 56);
#endif
      }
      return val;
    }

  }

  
  template <
    unsigned CHUNK_SIZE = 1024 * 1024,    // per thread
    unsigned CHUNKS = 16                  // in flight
  >
  class writer_t {
  public:
    static constexpr char gz_header[] = "\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x03";
    static constexpr unsigned chunks = CHUNKS;
    static constexpr unsigned chunk_size = CHUNK_SIZE;

    unsigned tot_chunk_size() const noexcept { return num_threads_ * chunk_size; }
    unsigned num_threads() const noexcept { return num_threads_; }
    
    // approximated bytes compressed, accurate after close    
    size_t tot_bytes() const noexcept { return tot_bytes_written_.load(std::memory_order_relaxed); }

    writer_t(const std::filesystem::path& output, std::shared_ptr<hahi::pool_t> pool, unsigned num_threads = -1) 
    : in_chunks_{chunks},
      num_threads_(num_threads),
      shpool_(pool)
    {
      if (num_threads_ == -1) num_threads_ = pool->num_threads();
      num_threads_ = std::clamp(num_threads_, 1u, pool->num_threads());
      auto gzout = std::ofstream(output, std::ios::binary);
      gzout.write(gz_header, sizeof(gz_header) - 1);
      closed_ = false;
      in_chunk_.reserve(tot_chunk_size());
      launch_compressor(std::move(gzout), shpool_, num_threads);
    }

    ~writer_t() {
      close(true);
      if (compressor_.joinable()) {
        compressor_.join();
      }
    }

    bool failed() const noexcept { return !!eptr_; }
    bool closed() const noexcept { return closed_; }
    
    void close(bool join = false) {
      if (closed_) return;
      closed_ = true;
      in_chunks_.emplace(std::move(in_chunk_)); // last, incomplete chunk
      if (join && compressor_.joinable()) {
        compressor_.join();
      }
    }

    // adds additional newline to the output
    template <typename T>
    void puts(const T& val) {
      if (closed_) throw std::runtime_error("fastq_writer: attempt to write into closed stream");
      if constexpr (std::is_convertible_v<T, std::string> ||
                    std::is_convertible_v<T, str_view>) {
        do_put<true>(val);
      }
      else {
        for (auto& str : val) do_put<true>(str);
      }
    }

    // doesn't add additional newline to the output
    template <typename T>
    void put(const T& val) {
      if (closed_) throw std::runtime_error("fastq_writer: attempt to write into closed stream");
      if constexpr (std::is_convertible_v<T, std::string> ||
                    std::is_convertible_v<T, str_view>) {
        do_put<false>(val);
      }
      else {
        for (auto& str : val) do_put<false>(str);
      }
    }

  private:
    template <bool newline>
    void do_put(str_view str) {
      assert((str.length() + newline) < tot_chunk_size());
      auto avail = tot_chunk_size() - in_chunk_.size();
      if (avail < str.length() + newline) {
        in_chunk_.insert(in_chunk_.end(), str.cbegin(), str.cbegin() + avail);
        assert(in_chunk_.length() == tot_chunk_size());
        in_chunks_.emplace(std::move(in_chunk_));
        str.remove_prefix(avail);
      }
      in_chunk_.insert(in_chunk_.cend(), str.cbegin(), str.cend());
      if constexpr (newline) {
        in_chunk_.push_back('\n');
      }
    }

    void launch_compressor(std::ofstream&& gzout, auto shpool, unsigned nun_threads) {
      // this is a rip-off of Mark Adler's pigz code: https://zlib.net/pigz/
      // 
      compressor_ = std::thread([&, gzout = std::move(gzout), pool = shpool_.get(), nt = num_threads_]() mutable {
        using cfuture = std::future<std::pair<char* /* out buf */, uint32_t /* avail */>>;
        auto cf = std::vector<cfuture>{};

        // set up numtreads z_streams
        const size_t gz_buffer = (4 * chunk_size) / 3;
        std::vector<z_stream> strms{size_t(nt)};
        for (auto& strm : strms) {
          std::memset(&strm, 0, sizeof(strm));
        }
        try {
          auto raw_out = std::unique_ptr<char[]>(new char[nt * gz_buffer + 4096]);
          auto out = (char*)(void*)((uintptr_t(raw_out.get()) + 4095) & ~4095);  // align to page size
          for (auto& strm : strms) {
            z_stream_init(strm);
          }
          uint32_t crc = crc32(0L, nullptr, 0);
          uint64_t tot_bytes = 0;
          auto prep_strm = [&](int idx, uint32_t avail_in, char* buf) -> unsigned {
            auto tmp = std::min(avail_in, chunk_size);
            auto& strm = strms[idx];
            z_stream_reset(strm);
            strm.avail_in = tmp;
            strm.next_in = (unsigned char*)buf + idx * chunk_size;
            strm.avail_out = gz_buffer;
            strm.next_out = (unsigned char*)out + idx * gz_buffer;
            strm.data_type = Z_SYNC_FLUSH;  // abuse, fixed later
            return avail_in - tmp;
          };
          for (bool last = false; !last;) {
            auto buf = in_chunks_.pop();
            tot_bytes += buf.size();
            tot_bytes_written_.store(tot_bytes, std::memory_order_relaxed);
            // allocate work to compress threads
            auto avail_in = static_cast<uint32_t>(buf.size());
            last = avail_in < tot_chunk_size();    // stop condition
            int nct = 0;  // used threads
            do  {
              avail_in = prep_strm(nct++, avail_in, buf.data());
            } while (avail_in);
            if (last) strms[nct - 1].data_type = Z_FINISH;  // abuse, fixed later
            cf.clear();
            for (auto j = 0; j < nct; ++j) {
              cf.emplace_back(
                pool->async([](z_stream* strm) {
                  const auto out = strm->next_out;
                  const int flush = std::exchange(strm->data_type, 2);  // fix abuse
                  deflate(strm, flush);
                  assert(strm->avail_in == 0);   // all input consumed
                  return std::pair<char*, uint32_t>{ (char*)out, strm->avail_out };
                }, &strms[j])
              );
            }
            crc = crc32(crc, (unsigned char*)buf.data(), buf.size());
            for (auto j = 0; j < nct; ++j) {
              // write deflate chunks in order
              auto [cout, avail] = cf[j].get();
              gzout.write((const char*)cout, gz_buffer - avail);
            }
            gzout.flush();
          }
          // write 8 byte gz footer (little endian)
          crc = htogz(crc);
          tot_bytes = htogz(tot_bytes);
          gzout.write((const char*)&crc, 4);
          gzout.write((const char*)&tot_bytes, 4);
          tot_bytes_written_.store(tot_bytes, std::memory_order_release);
        }
        catch (...) {
          eptr_ = std::current_exception();
        }
        for (auto& strm : strms) (void)deflateEnd(&strm);
        gzout.close();
      });
    }

    template <typename T> 
    using queue_t = hahi::concurrent_queue<T>;
    using buffer_t = std::string;

    buffer_t in_chunk_;              // current input buffer used by put functions
    queue_t<buffer_t> in_chunks_;    // populated by put functions, consumed by compress thread
    std::exception_ptr eptr_;
    bool closed_ = true;
    unsigned num_threads_ = 0;
    std::shared_ptr<hahi::pool_t> shpool_;
    std::atomic<size_t> tot_bytes_written_ = 0;
    std::thread compressor_;
  };

}
