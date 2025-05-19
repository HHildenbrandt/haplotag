/*
 * Copyright (c) 2023 Hanno Hildenbrandt
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
 *
 */

#ifndef HAHI_MUTEX_HPP_INCLUDED
#define HAHI_MUTEX_HPP_INCLUDED
#pragma once

#include <atomic>
#include <mutex>

#if defined(_MSC_VER)
  #include <intrin.h>
  #if defined(_M_AMD64) || defined(_M_IX86)
    #pragma intrinsic(_mm_pause)
    #define hahi_spin_pause()  _mm_pause()
  #endif
#elif defined(__x86_64__) || defined(__i386__)
  #define hahi_spin_pause()  __asm__ __volatile__ ("pause")
#elif defined(__arm__)
  #ifdef __CC_ARM
    #define hahi_spin_pause()  __yield()
  #else
    #define hahi_spin_pause()  __asm__ __volatile__ ("yield")
  #endif
#endif


namespace hahi {


  inline void spin_pause() noexcept {
    hahi_spin_pause();
  }


  // lightweight, unfair spin-mutex
  // good performance under very low congestion
  // high power consumption under congestion due to busy-loop
  class spin_lock
  {
    mutable std::atomic_flag flag_ = ATOMIC_FLAG_INIT;

  public:
    void lock() const noexcept {
      while (flag_.test_and_set(std::memory_order_acquire)) {
        spin_pause();
      }
    }

    [[nodiscard]]
    bool try_lock() const noexcept {
      return !flag_.test_and_set(std::memory_order_acquire);
    }

    void unlock() const noexcept {
      flag_.clear(std::memory_order_release);
    }
  };


  // lightweight, unfair spin-mutex
  // good performance under very low congestion
  // may back-of to yield and sleep
  class spin_mutex {
  public:
    void lock() const noexcept {
      while (flag_.test_and_set(std::memory_order_acquire)) {
        // wait for flag_ to change to false
        flag_.wait(true, std::memory_order_relaxed);
      }
    }

    [[nodiscard]]
    bool try_lock() const noexcept {
      return !flag_.test_and_set(std::memory_order_acquire);
    }

    void unlock() const noexcept {
      flag_.clear(std::memory_order_release);
      flag_.notify_one();
    }

  private:
    mutable std::atomic_flag flag_ = ATOMIC_FLAG_INIT;
  };


  // lightweight, fair spin-mutex
  // poor performance under high congestion
  // may back-of to yield and sleep
  // no support for try_lock
  class ticket_mutex {
  public:
    void lock() const noexcept {
      const auto ticket = ticket_.fetch_add(1, std::memory_order_acquire);
      while (true) {
        const auto served = served_.load(std::memory_order_acquire);
        if (ticket == served) return;
        // wait for served_ to change
        served_.wait(served, std::memory_order_relaxed);
      }
    }

    void unlock() const noexcept {
      served_.fetch_add(1, std::memory_order_release);
      served_.notify_all();
    }

  private:
    mutable std::atomic<int> ticket_ = ATOMIC_VAR_INIT(0);
    mutable std::atomic<int> served_ = ATOMIC_VAR_INIT(0);
  };

}

#endif
