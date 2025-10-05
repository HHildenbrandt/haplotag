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

#ifndef HAHI_POOL_HPP_INCLUDED
#define HAHI_POOL_HPP_INCLUDED
#pragma once
 
#include <bit>
#include <vector>
#include <memory>
#include <algorithm>     // std::clamp
#include <bitset>
#include "device.hpp"


namespace hahi {

  // models a simple fixed-sized thread pool.
  //
  // provides limited concurrent forward progress guarantees!
  // avoid infinite tasks as they will dead-lock the whole pool eventually!
  //
  class pool_t {
  public:
    // arbitrary limit in multiplies of 64
    static constexpr unsigned max_threads = 256;

    // creates pool with num_threads clamped to a hardware concurrency
    explicit pool_t(unsigned num_threads = -1) 
    : sem_{std::clamp(num_threads, 1u, std::thread::hardware_concurrency())} {
      num_threads = std::clamp(num_threads, 1u, std::thread::hardware_concurrency());
      if (num_threads > max_threads) {
        throw std::runtime_error("Number of threads exceeds implementation limit");
      }
      auto fit = free_list_;
      auto f64 = num_threads;
      for (; f64 >= 64; f64 -= 64, ++fit) *fit = -1;
      if (f64) *fit = (1ull << f64) - 1;
      for (unsigned i = 0; i < num_threads; ++i) {
        devices_.emplace_back(new device_t{1 + 1});   // 1 work + 1 release task
      }
    }

    unsigned num_threads() const noexcept { return devices_.size(); }

    // returns number of available jobs
    int avail() const noexcept { 
      std::lock_guard<std::mutex> _{mutex_};
      int bits = 0;
      for (size_t i = 0; i < sizeof(free_list_); ++i) {
        bits += std::popcount(free_list_[i]);
      }
      return bits; 
    }
  
    // returns number of running jobs
    int busy() const noexcept { return num_threads() - avail(); }

    // submits job to pool.
    // returns std::future
    template <typename Fun, typename... Args>
    auto async(Fun&& fun, Args&&...args) const {
      sem_.acquire();   // wait for idle devices
      int f64 = 0;
      int bit = 0;
      {
        std::lock_guard<std::mutex> _(mutex_);
        while (0 == free_list_[f64]) ++f64;
        bit = std::countr_zero(free_list_[f64]);
        free_list_[f64] &= ~(1ull << bit);
      }
      auto future = devices_[bit]->enqueue(std::forward<Fun>(fun), std::forward<Args>(args)...);
      devices_[bit]->enqueue_detach([&, bit = bit, f64 = f64]() noexcept { 
        std::lock_guard<std::mutex> _(mutex_);
        free_list_[f64] |= (1ull << bit);
        sem_.release(1);
      });
      return future;
    }

  private:
    mutable std::counting_semaphore<> sem_;
    mutable std::mutex mutex_;
    mutable uint64_t free_list_[max_threads >> 6];    // bitset
    std::vector<std::unique_ptr<device_t>> devices_;
  };

}

#endif
