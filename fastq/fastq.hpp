/* 
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

#include <string_view>
#include <algorithm>
#include <memory>


namespace fastq {

  using str_view = std::string_view;
  using chunk_ptr = std::shared_ptr<char[]>;


  // returns exhaustiv str_view if [pos, pos+length()) is out of bounds
  inline str_view max_substr(str_view str, 
                             size_t pos) noexcept {
    auto first = str.data() + pos;
    auto last = str.data() + str.length();
    return { std::min(first, last), last };
  }


  // returns exhaustiv str_view if [pos, pos+count) is out of bounds
  inline str_view max_substr(str_view str, 
                             size_t pos,
                             size_t count) noexcept {
    auto first = str.data() + pos;
    auto last = str.data() + str.length();
    return { std::min(first, last), std::min(first + count, last) };
  }

}
