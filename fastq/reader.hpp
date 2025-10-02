/* fastq/reader.hpp
 *
 * Copyright (c) 2025 Hanno Hildenbrandt <h.hildenbrandt@rug.nl>
 */

#pragma once

#include <cassert>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <stdexcept>
#include <limits>
#include <memory>
#include <filesystem>
#include <thread>
#include <future>
#include <atomic>
#include <utility>
#include <array>
#include <list>
#include <span>
#include <string_view>
#include <device/queue.hpp>
#include <device/mutex.hpp>
#include "fastq.hpp"


namespace fastq {


  // blob handed out by reader_t<...>::operator()
  //
  // buf ->          |..............| }
  //                       ...        |- undefined spare
  //                 |..............| }  
  // buf + window -> |..............|                }
  //                 |..............|                |
  //                       ...                       |- cv()
  //                       ...                       |
  //                 |..............| <- buf + size  }
  //                       ...
  //                 |..............| <- buf + reader_t<...>::chunk_size
  // 
  struct chunk_t {
    chunk_ptr buf;
    
    size_t size = 0;      // available characters
    size_t window = 0;    // offset first character
    bool last = false;    // last chunk available from reader

    char* data() const noexcept { return buf.get() + window; }
    str_view cv() const noexcept { return str_view{data(), size }; }
  };

  
  namespace detail {

    // asynchronous wrapper around `zlib::gzread`
    // as such, accepts uncompressed files too
    template <
      typename Allocator,
      size_t Window = 16 * 1024,          // shall be bigger than max item size
      size_t ChunkSize = 1024 * 1024,     // max. chunk size (exclusive window. padding)
      unsigned Chunks = 16,               // queue depth, chunks in flight
      unsigned GzBuffer = 128 * 1024
    >
    class reader_t {
    public:
      static_assert(Window < (ChunkSize >> 4));
      static_assert(ChunkSize < std::numeric_limits<int>::max());   // zlib limitation
      static constexpr size_t window = Window;
      static constexpr size_t chunk_size = ChunkSize;
      static constexpr unsigned chunks = Chunks;
      static constexpr unsigned gz_buffer = GzBuffer;
      using allocator_t = Allocator;

      reader_t() = default;
      reader_t(reader_t&&) = default;
      reader_t& operator=(reader_t&&) = default;
      
      explicit reader_t(const std::filesystem::path& path) : path_(path) {
        if (nullptr == (gzin_ = zng_gzopen(path.string().c_str(), "rb"))) {
          throw std::runtime_error(std::string("fastq::reader_t: failed to open input file \'") + path.string() + '\'');
        }
        zng_gzbuffer(gzin_, gz_buffer);
        deflate_ = std::jthread([&, gzin = gzin_](std::stop_token stok) {
          try {
            while (!stok.stop_requested()) {
              auto buf =  chunk_ptr(static_cast<char*>(alloc_.alloc(chunk_size + window)), &allocator_t::free);
              const auto avail = static_cast<size_t>(zng_gzread(gzin, buf.get() + window, static_cast<unsigned>(chunk_size)));
              if (avail == size_t(-1)) {  // error
                throw -1;
              }
              const bool last = avail < chunk_size;
              chunks_.push(chunk_t{ .buf = std::move(buf), .size = avail, .window = window, .last = last});
              if (last) {   // eof
                break;
              }
            }
          }
          catch (...) { // sink exception
            fail_.store(true, std::memory_order_release);
          }
          chunks_.emplace(nullptr, 0);    // sentinel
        });
      }

    public:
      ~reader_t() {
        // gracefully end worker threads if necessary
        deflate_.request_stop();
        if (deflate_.joinable()) {
          while (chunks_.try_pop().has_value()) ;   // deplete file queue. allow reader_ to push sentinel
          deflate_.join();
        }
        zng_gzclose(gzin_);
      }

      // bytes deflated
      size_t tot_bytes() const noexcept { return tot_bytes_; }
      bool failed() const noexcept { return fail_.load(std::memory_order_acquire); }
      bool eof() const noexcept { return eof_; }
      const std::filesystem::path& path() const noexcept { return path_; }
      allocator_t& allocator() { return alloc_; }

      // returns new chunk or an empty chunk_t if eof() == true
      chunk_t operator()() {
        if (!eof_) {
          auto chunk = chunks_.pop(); 
          tot_bytes_ += chunk.size;
          eof_ = chunk.last | fail_.load(std::memory_order_acquire);
          return chunk;
        }
        return {};
      }

    private:
      mutable hahi::concurrent_queue<chunk_t> chunks_{chunks};
      mutable std::atomic<bool> fail_{false};
      size_t tot_bytes_ = 0;
      bool eof_ = false;
      allocator_t alloc_;
      std::jthread deflate_;
      gzFile gzin_ = nullptr;
      const std::filesystem::path path_;
    };


    struct default_allocator {
      static void* alloc(size_t bytes) { 
        void* ptr = std::malloc(bytes);
        if (nullptr == ptr) throw std::bad_alloc{}; 
        return ptr;
      }
        
      static void free(void* ptr) noexcept { std::free(ptr); } 
    };

  }


  // asynchronous wrapper around `zlib::gzread`
  using reader_t = detail::reader_t<detail::default_allocator>;

}
