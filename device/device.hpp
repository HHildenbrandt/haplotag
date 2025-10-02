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

#ifndef HAHI_DEVICE_HPP_INCLUDED
#define HAHI_DEVICE_HPP_INCLUDED
#pragma once

#include <type_traits>
#include <deque>
#include <mutex>
#include <future>
#include <thread>
#include <stop_token>
#include <semaphore>
#include "mutex.hpp"    // spin_mutex
#include "queue.hpp"
#include <functional>
#if !__cpp_lib_move_only_function 
# if __has_include(<function2/function2.hpp>)
  // https://github.com/Naios/function2.git
# include <function2/function2.hpp>
# endif
#endif


namespace hahi {

  namespace detail {

#if __cpp_lib_move_only_function
    using task_function = std::move_only_function<void()&&>;
#elif defined(FU2_INCLUDED_FUNCTION2_HPP_)
    // https://github.com/Naios/function2.git
    using task_function = fu2::unique_function<void()&&>;
#else
    using task_function = std::function<void()>;
#endif


    template <typename T>
    inline std::decay_t<T> decay_copy(T&& v) { 
      return std::forward<T>(v); 
    }

    // wraps Fun(Args...) into a nullary one-shot lambda. Lambda captures by value.
    template <typename... Args, typename Fun> 
    inline auto wrap(Fun&& fun, Args &&...args) {
      if constexpr (sizeof...(Args) == 0) {
        return fun;
      }
      else {
        return [fun = decay_copy(std::forward<Fun>(fun)), ... args = decay_copy(std::forward<Args>(args))]() mutable -> decltype(auto) {
          return std::invoke(std::move(fun), std::move(args)...);
        };
      }
    }


    // resembles std::packaged_task<R() &&>
    template <std::invocable Fun> 
    class one_shot_task  {
    public:
      // Stores a copy of the function
      one_shot_task(Fun fun) : fun_(std::move(fun)) {}

      std::future<std::invoke_result_t<Fun>> get_future() { return promise_.get_future(); }

      // invocable via std::invoke(std::move(one_shot_task))
      void operator()() && {
        try {
          if constexpr (!std::is_same_v<void, std::invoke_result_t<Fun>>) {
            promise_.set_value(std::invoke(std::move(fun_)));
          }
          else {
            std::invoke(std::move(fun_));
            promise_.set_value();
          }
        }
        catch (...) {
          promise_.set_exception(std::current_exception());
        }
      }

    private:
      std::promise<std::invoke_result_t<Fun>> promise_;
      Fun fun_;
    };

  }  // namespace detail


  // models a single-threaded thread-pool
  class device_t {
    using queue_t = concurrent_queue<detail::task_function>;

  public:
    explicit device_t(unsigned max_pending) : queue_{max_pending} {
      thread_ = std::jthread([&](std::stop_token stoken) {
        do {
          std::invoke(queue_.template pop<typename queue_t::explicit_release>());
          queue_.release();   // signal work completion
        } while (!stoken.stop_requested());
      });
    };

    ~device_t() {
      enqueue([](){}).get();    // send & wait on 1st sentinel (flush queue)
      thread_.request_stop();
      enqueue_detach([](){});   // send 2nd sentinel (fool task_sem_)
    };

    // enque job
    // returns the job's std::future
    template <typename Fun, typename... Args>
    [[nodiscard]] std::future<std::invoke_result_t<std::decay_t<Fun>, std::decay_t<Args>...>> enqueue(Fun&& fun, Args&&...args) {
      auto task = detail::one_shot_task(detail::wrap(std::forward<Fun>(fun), std::forward<Args>(args)...));
      auto future = task.get_future();
      queue_.push(std::move(task));
      return future;
    }

    // Doesn't return a future
    // Doesn't handle exceptions. If Fun throws, std::terminate is called
    template <typename Fun, typename... Args>
    void enqueue_detach(Fun&& fun, Args&&...args) {
      static_assert(std::is_same_v<void, std::invoke_result_t<std::decay_t<Fun>, std::decay_t<Args>...>>, "detached function must return void.");
      queue_.push(detail::wrap(std::forward<Fun>(fun), std::forward<Args>(args)...));
    }

  private:
    queue_t queue_;
    std::jthread thread_;
  };

}

#endif
