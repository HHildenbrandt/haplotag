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
#include "device.hpp"


namespace hahi {

  // models a simple thread pool.
  //
  // provides limited concurrent forward progress guarantees!
  // avoid infinite tasks as they will dead-lock the whole pool eventually!
  //
  class pool_t {
  public:
    static constexpr unsigned max_threads = 256;   // current 'limitation'

    explicit pool_t(unsigned num_threads = -1) 
    : sem_{std::clamp(num_threads, 1u, std::thread::hardware_concurrency())} {
      num_threads = std::clamp(num_threads, 1u, std::thread::hardware_concurrency());
      if (num_threads > 256) {
        throw std::runtime_error("Dude, what are you running? More than 256 threads are not supported yet");
      }
      const auto n64 = num_threads >> 6;
      for (auto f64 = 0; f64 < n64; ++f64) free_list_[f64] = -1;
      free_list_[n64] = ((num_threads - (n64 << 6)) == 64) ? uint64_t(-1) : ((1ull << (num_threads - (n64 << 6))) - 1);
      for (unsigned i = 0; i < num_threads; ++i) {
        devices_.emplace_back(new device_t{1 + 1});   // 1 work + 1 release task
      }
    }

    unsigned num_threads() const noexcept { return devices_.size(); }
    
    int busy() const noexcept { 
      std::lock_guard<std::mutex> _{mutex_};
      int f64 = 0; while (0 == free_list_[f64]) ++f64;
      return num_threads() - (f64 * 64) + std::popcount(free_list_[f64]); 
    }
  
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
    mutable uint64_t free_list_[4];    // bitset
    std::vector<std::unique_ptr<device_t>> devices_;
  };

}

#endif
