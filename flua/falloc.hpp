// falloc.h allocator for flua

// Copyright (C) 2024 Hanno Hildenbrandt
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files(the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and /or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions :
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef FLUA_FALLOC_HPP_INCLIDED
#define FLUA_FALLOC_HPP_INCLIDED
#pragma once


#include <cassert>
#include <cstdlib>
#include <new>        // std::bad_alloc
#include <stdexcept>
#include <limits>
#include <memory>
#include <bit>
#include <array>
#include <utility>
#include <cmath>
#include <cstring>
#include <cstddef>

#if defined(FLUA_USE_WINDOWS)
# include <windows.h>
# include <sysinfoapi.h>
# include <memoryapi.h>
#else /* Linux, posix, darwin */
# include <sys/mman.h>
# include <unistd.h>
#endif


namespace flua {


  struct l_alloc {
    static void* frealloc(void* /* ud */, void *ptr, size_t /* osize */, size_t nsize) {
      if (nsize == 0) {
        std::free(ptr);
        return nullptr;
      }
      else {
        return std::realloc(ptr, nsize);
      }
    }
  };


  namespace detail {

    template <size_t align = alignof(std::max_align_t)>
    constexpr size_t aligned_bytes(size_t bytes) noexcept {
      static_assert(std::has_single_bit(align), "'align' shall be power of two");
      return (bytes + align - 1) & ~(align - 1);
    }

 
#if defined(FLUA_USE_WINDOWS)

    inline void* virtual_alloc(size_t nbytes) noexcept {
      return VirtualAlloc(nullptr, nbytes, MEM_RESERVE, PAGE_READWRITE);
    }


    inline void* virtual_commit(void* start, size_t nbytes) noexcept {
      return VirtualAlloc(start, nbytes, MEM_COMMIT, PAGE_READWRITE);
    }


    inline void virtual_decommit(void* ptr, size_t nbytes) noexcept {
      VirtualFree(ptr, nbytes, MEM_DECOMMIT);
    }


    inline void virtual_free(void* ptr, size_t /* nbytes */) noexcept {
      VirtualFree(ptr, 0, MEM_RESET);
    }


    inline size_t page_size() {
      SYSTEM_INFO si;
      GetSystemInfo(&si);
      return static_cast<size_t>(si.dwPageSize);
    }

#else /* Linux | POSIX */

    inline void* virtual_alloc(size_t nbytes) noexcept {
      auto ptr = mmap(NULL, nbytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
      return (MAP_FAILED == ptr) ? nullptr : ptr;
    }


    inline void* virtual_commit(void* /*start*/, size_t /*nbytes*/) noexcept {
      // no-op in Linux
      return reinterpret_cast<void*>(1);
    }


    // `ptr` must be aligned to page_size
    inline void virtual_decommit(void* ptr, size_t nbytes) noexcept {
#if defined(FLUA_USE_LINUX)      
      madvise(ptr, nbytes, MADV_DONTNEED);
#else /* POSIX */
      posix_madvise(ptr, nbytes, POSIX_MADV_DONTNEED);
#endif     
    }


    inline void virtual_free(void* ptr, size_t nbytes) noexcept {
      munmap(ptr, nbytes);
    }

    
    inline size_t page_size() {
      return getpagesize();
    }

#endif

  } // namespace detail


  // Reserves `MaxMB` MB in virtual address space and handles
  // out memory chunks of size `BlkSize` aligned to `BlkSize`
  template <
    size_t MaxMB,
    size_t BlkSize,
    bool Decommit
  >
  class vchunk_allocator_t {
  public:
    static constexpr size_t ReqBytes = 1024 * 1024 * MaxMB + BlkSize;
    static constexpr size_t ChunkSize = 64 * BlkSize;
    static constexpr size_t NumChunks = std::max<size_t>(ReqBytes / ChunkSize, 1);
    static constexpr size_t ReservedBytes = (NumChunks + 1) * ChunkSize;

    vchunk_allocator_t(const vchunk_allocator_t&) = delete;
    vchunk_allocator_t& operator=(const vchunk_allocator_t&) = delete;

    vchunk_allocator_t& operator=(vchunk_allocator_t&& rhs) noexcept {
      bitset_ = rhs.bitset_; rhs.bitset_.fill(-1);
      chunk_ = std::exchange(rhs.chunk_, 0);
      first_ = std::exchange(rhs.first_, nullptr);
      raw_ = std::exchange(rhs.raw_, nullptr);
      return *this;
    }

    vchunk_allocator_t(vchunk_allocator_t&& rhs) noexcept {
      *this = std::move(rhs);
    }

    vchunk_allocator_t() {
      if constexpr (Decommit) {
        if (0 != (BlkSize % detail::page_size())) {
          throw std::runtime_error("vchunk_allocator_t: `BlkSize` shall be a multiple of page size");
        }
      }
      raw_ = detail::virtual_alloc(ReservedBytes);
      if (nullptr == raw_) {
        throw std::bad_alloc();
      }
      bitset_.fill(raw_ ? 0 : -1);
      first_ = (void*)((uintptr_t(raw_) + BlkSize - 1) & ~(BlkSize - 1));
    }

    ~vchunk_allocator_t() noexcept {
      if (raw_) {
        // Debug hint: expect a single bit set in bitset_
        detail::virtual_free(raw_, ReservedBytes);
      }
    }

    void* alloc() noexcept {
      auto chunk = chunk_;
      auto bits = bitset_[chunk];
      for (; -1 == bits; bits = bitset_[++chunk])
        ;
      if (chunk == bitset_.size() - 1) [[unlikely]] {
        return nullptr;     // out of memory, triggers LUA's emergency GC
      }
      const auto bit = std::countr_one(bits);
      const auto cptr = chunk_ptr(first_, chunk);
      void* const ptr = (void*)(cptr + BlkSize * bit);
      if (0 == bits) [[unlikely]] {
        if (nullptr == detail::virtual_commit((void*)cptr, ChunkSize)) [[unlikely]] {
          return nullptr;
        }
      }
      bitset_[chunk] |= (uint64_t(1) << bit);
      chunk_ = chunk;
      return ptr;
    }

    void free(void* blk) noexcept {
      const size_t blk_idx = (uintptr_t(blk) - uintptr_t(first_)) / BlkSize;
      const size_t chunk = blk_idx / 64;
      const auto bits = (bitset_[chunk] &= ~(uint64_t(1) << (blk_idx - 64 * chunk)));
      if constexpr (Decommit) {
        if (0 == bits) [[unlikely]] {
          detail::virtual_decommit((void*)chunk_ptr(first_, chunk), ChunkSize);
        }
      }
      chunk_ = std::min(chunk_, chunk);
    }

  private:
    static constexpr std::byte* chunk_ptr(void* first, size_t chunk) {
      return std::bit_cast<std::byte*>(first) + chunk * ChunkSize;
    }

    std::array<uint64_t, NumChunks + 1> bitset_;   // in-use list, +1: guard
    size_t chunk_ = 0;        // current allocation chunk
    void* first_ = nullptr;   // BlkSize aligned
    void* raw_ = nullptr;
  };


  // A memory allocator that works pretty well with the generational GC of Lua 5.4.
  // Consist of a swarm of stack-allocators of size BlkSize backed by vblock_allocator.
  template <
    size_t MaxMB = 1024,
    size_t BlkSize = 64 * 1024,
    bool Decommit = true,
    typename Fallback = flua::l_alloc
  >
  class allocator_t {
    struct alignas(std::max_align_t) block_t {
      uint32_t alloc;   // offset next allocation
      uint32_t n;       // number of allocations
    };

    static constexpr uint32_t block_capacity = static_cast<uint32_t>(BlkSize - sizeof(block_t));

    static void* blk_allocate(block_t* blk, uint32_t nbytes) noexcept {
      const auto ofs = blk->alloc;
      blk->alloc += nbytes;
      ++blk->n;
      return (void*)(std::bit_cast<std::byte*>(blk) + ofs);
    }

  public:
    allocator_t(const allocator_t&) = delete;
    allocator_t(allocator_t&&) = default;

    allocator_t() {
      if (nullptr == (head_ = alloc_block())) {
        throw std::bad_alloc{};
      }
    }

    /*
    ** From lmem.c:
    **
    ** About the realloc function:
    ** void *frealloc (void *ud, void *ptr, size_t osize, size_t nsize);
    ** ('osize' is the old size, 'nsize' is the new size)
    **
    ** - frealloc(ud, p, x, 0) frees the block 'p' and returns NULL.
    ** Particularly, frealloc(ud, NULL, 0, 0) does nothing,
    ** which is equivalent to free(NULL) in ISO C.
    **
    ** - frealloc(ud, NULL, x, s) creates a new block of size 's'
    ** (no matter 'x'). Returns NULL if it cannot create the new block.
    **
    ** - otherwise, frealloc(ud, b, x, y) reallocates the block 'b' from
    ** size 'x' to size 'y'. Returns NULL if it cannot reallocate the
    ** block to the new size.
    */
    static void* frealloc(void* ud, void* optr, size_t osize, size_t nsize) noexcept {
      auto self = static_cast<allocator_t*>(ud);
      osize = detail::aligned_bytes(osize);
      nsize = detail::aligned_bytes(nsize);
      if (nsize) {
        return (nullptr == optr)
          ? self->allocate(nsize)
          : self->reallocate(optr, osize, nsize);
      }
      if (optr) [[likely]] {
        self->deallocate(optr, osize);
      }
      return nullptr;
    }

    static void* alloc(void* ud, size_t bytes) noexcept {
      bytes = detail::aligned_bytes(bytes);
      return (bytes) ? static_cast<allocator_t*>(ud)->allocate(bytes) : nullptr;
    }

    static void free(void* ud, void* ptr, size_t bytes) noexcept {
      if (ptr) [[likely]] {
        static_cast<allocator_t*>(ud)->deallocate(ptr, bytes);
      }
    }

    
  private:
    block_t* alloc_block() noexcept {
      if (void* ptr = vchunks_.alloc()) {
        block_t* blk = new (ptr) block_t{};
        blk->alloc = detail::aligned_bytes(sizeof(block_t));
        blk->n = 0;
        return blk;
      }
      return nullptr;
    }

    void free_block(block_t* blk) noexcept {
      assert(0 == blk->n);
      blk->alloc = detail::aligned_bytes(sizeof(block_t));   // reclaim storage
      if (blk != head_) [[likely]] {
        // not head_, free block
        blk->~block_t();  // trivial dtor
        vchunks_.free(blk);
      }
    }

    // blk(ptr)
    constexpr block_t* lookup(void* ptr) const noexcept {
      return std::bit_cast<block_t*>(uintptr_t(ptr) & ~(BlkSize - 1));
    }

    void* allocate(uint32_t nbytes) noexcept {
      if (nbytes > block_capacity) [[unlikely]] {
        return Fallback::frealloc(&fallback_, nullptr, 0, nbytes);
      }
      if ((head_->alloc + nbytes) > block_capacity) [[unlikely]] {
        auto head = alloc_block();
        if (nullptr == head) [[unlikely]] {
          return nullptr;
        }
        head_ = head;
      }
      return blk_allocate(head_, nbytes);
    }

    void deallocate(void* ptr, uint32_t ubytes) noexcept {
      if (ubytes > block_capacity) [[unlikely]] {
        Fallback::frealloc(&fallback_, ptr, ubytes, 0);
        return;
      }
      block_t* const blk = lookup(ptr);
      if (ptr == (void*)(std::bit_cast<std::byte*>(blk) + (blk->alloc - ubytes))) {
        blk->alloc -= ubytes; // reclaim
      }
      if (0 == --blk->n) [[unlikely]] {
        free_block(blk);
      }
    }

    void* reallocate(void* ptr, uint32_t obytes, uint32_t nbytes) noexcept {
      if (obytes > block_capacity) [[unlikely]] {
        if (nbytes > block_capacity) [[likely]] {
          return Fallback::frealloc(&fallback_, ptr, obytes, nbytes);
        }
        // take ownership of shrunken memory block
        void* nptr = allocate(nbytes);
        if (nptr) [[likely]] {
          std::memcpy(nptr, ptr, std::min(obytes, nbytes));
          Fallback::frealloc(&fallback_, ptr, obytes, 0);
        }
        return nptr;
      }
      block_t* const blk = lookup(ptr);
      if (ptr == (void*)(std::bit_cast<std::byte*>(blk) + (blk->alloc - obytes))) {
        blk->alloc -= obytes; // reclaim
        if ((blk->alloc + nbytes) < block_capacity) [[likely]] {
          blk->alloc += nbytes;
          return ptr;
        }
      }
      if (obytes >= nbytes) [[unlikely]] {
        return ptr;
      }
      if (void* const nptr = allocate(nbytes); nptr) [[likely]] {
        std::memcpy(nptr, ptr, std::min(obytes, nbytes));
        if (0 == --blk->n) [[unlikely]] {
          free_block(blk);
        }
        return nptr;
      }
      return nullptr;
    }

    block_t* head_ = nullptr;    // current allocation node
    vchunk_allocator_t<MaxMB, BlkSize, Decommit> vchunks_;
    Fallback fallback_;
  };


} // namepsace flua


#endif
