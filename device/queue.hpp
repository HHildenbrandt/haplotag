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

#ifndef HAHI_QUEUE_HPP_INCLUDED
#define HAHI_QUEUE_HPP_INCLUDED
#pragma once

#include <mutex>
#include <semaphore>
#include <optional>
#include <memory>
#include <queue>


namespace hahi {


  // A fixed-size concurrent queue
  template <typename T, typename Mutex = std::mutex>
  class concurrent_queue
  {
  public:
    struct implicit_release {};
    struct explicit_release {};
    using value_type = T;

    explicit concurrent_queue(size_t max_size) 
    : in_sem_(max_size),
      out_sem_(0),
      queue_(new T[max_size]),
      max_size_(max_size)
    {}
    
    ~concurrent_queue() {}

    size_t max_size() const noexcept { return max_size_; }

    template <typename... Args>
    void emplace(Args&&... args) const {
      // syntactic sugar - performing move-assignment
      push(value_type(std::forward<Args>(args)...));
    }

    void push(value_type&& val) const {
      in_sem_.acquire();    // wait if queue is full
      {
        std::lock_guard<Mutex> _(mutex_);
        auto back = back_++ % max_size_;
        queue_[back] = std::move(val);
      }
      out_sem_.release(1);   // signal items available
    }

    bool try_push(value_type& val) const {
      if (in_sem_.try_acquire()) {
        {
          std::lock_guard<Mutex> _(mutex_);
          auto back = back_++ % max_size_;
          queue_[back] = std::move(val);
        }
        out_sem_.release(1);   // signal items available
        return true;
      }
      return false;
    }

    // ReleasePolicy == implicit_release: releases queue
    // ReleasePolicy == explicit_release: user must call releases()
    template <typename ReleasePolicy = implicit_release>
    void get(value_type& val) const {
      out_sem_.acquire();   // wait for items;
      val = dequeue<ReleasePolicy>();
    }

    // ReleasePolicy == implicit_release: releases queue
    // ReleasePolicy == explicit_release: user must call releases()
    template <typename ReleasePolicy = implicit_release>
    value_type pop() const {
      out_sem_.acquire();   // wait for items;
      return dequeue<ReleasePolicy>();
    }

    // ReleasePolicy == implicit_release: releases queue
    // ReleasePolicy == explicit_release: user must call releases()
    template <typename ReleasePolicy = implicit_release>
    std::optional<value_type> try_pop() const {
      if (out_sem_.try_acquire()) {
        return dequeue<ReleasePolicy>();
      }
      return {};
    }

    // explicit release
    // must be called after pop<explicit_release>() or if
    // try_pop<explicit_relaese> returns value
    void release() const {
      in_sem_.release(1);    // signal item consumption
    }

    bool try_acquire() const noexcept { return in_sem_.try_acquire(); }

  private:
    template <typename ReleasePolicy>
    T dequeue() const {
      std::unique_lock<Mutex> lock(mutex_);
      auto front = front_++ % max_size_;
      auto val = std::move(queue_[front]);
      lock.unlock();
      if constexpr (std::is_same_v<ReleasePolicy, implicit_release>) {
        in_sem_.release(1);    // signal item consumption
      }
      return val;
    }

    mutable std::counting_semaphore<> in_sem_;
    mutable std::counting_semaphore<> out_sem_;
    mutable Mutex mutex_;     // protects the queue
    mutable size_t front_ = 0;
    mutable size_t back_ = 0;
    mutable std::unique_ptr<T[]> queue_;
    const size_t max_size_;
  };

}

#endif // HAHI_QUEUE_HPP_INCLUDED
