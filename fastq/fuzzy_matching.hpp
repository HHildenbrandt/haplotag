/* fastq/fuzzy_matching.hpp
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
#include <string_view>
#include <numeric>  // iota, max
#include <algorithm>
#include <unordered_map>
#include "fastq.hpp"
#include "barcode.hpp"


namespace fastq {


  // generic Levensthein edit distance
  inline int edit_distance(str_view a, str_view b) {
    if (a.size() > b.size()) std::swap(a, b);
    // remove common prefix. removing suffix not worth it performance wise...
    while (!a.empty() && (a[0] == b[0])) { a.remove_prefix(1), b.remove_prefix(1); }
    int D[16 * 1024];  // just big enough
    const auto m = static_cast<int>(std::size(a));
    const auto n = static_cast<int>(std::size(b));
    std::iota(D, D + m + 1, 0);
    for (auto i = 1; i <= n; ++i) {
      auto tmp = std::exchange(D[0], i);  // tmp <- L(i-1,0), L(i,0) <- i
      for (auto j = 1; j <= m; ++j) {
        if (a[i - 1] != b[j - 1]) {
          tmp = std::min(D[j], std::min(D[j - 1], tmp)) + 1;
        }
        std::swap(tmp, D[j]);   // d[j] <- L(i,j)
      }
    }
    return D[m];
  }
  

  // bounded Levensthein edit distance, early exit if ed > min_ed
  inline int edit_distance(str_view a, str_view b, int min_ed) {
    if (a.size() > b.size()) std::swap(a, b);
    // remove common prefix. removing suffix not worth it performance wise...
    while (!a.empty() && (a[0] == b[0])) { a.remove_prefix(1), b.remove_prefix(1); }
    int D[16 * 1024];  // just big enough
    const auto m = static_cast<int>(std::size(a));
    const auto n = static_cast<int>(std::size(b));
    std::iota(D, D + m + 1, 0);
    for (auto i = 1; i <= n; ++i) {
      auto dmin = D[0];
      auto tmp = std::exchange(D[0], i);  // tmp <- L(i-1,0), L(i,0) <- i
      for (auto j = 1; j <= m; ++j) {
        if (a[i - 1] != b[j - 1]) {
          tmp = std::min(D[j], std::min(D[j - 1], tmp)) + 1;
        }
        dmin = std::min(dmin, tmp);
        std::swap(tmp, D[j]);   // d[j] <- L(i,j)
      }
      if (dmin > min_ed) {
        return dmin;    // can only get worse...
      }
    }
    return D[m];
  }
  

  enum ReadType{
    unclear,        // multiple occurrences of same min ed
    correct,        // exact match, ed == 0
    corrected,      // unique min ed
    max_read_type
  };


  struct match_t {
    const barcode_t& bc;
    int idx = -1;         // index into bc, undefined if false == valid()
    int ed = -1;          // edit distance, undefined if false == valid()
    ReadType read_type = ReadType::unclear;
    
    str_view tag() const { 
      return (read_type != unclear) ? bc[idx].tag() : bc.unclear_tag(); 
    }
    
    str_view code() const {
      return (read_type != unclear) ? bc[idx].code() : str_view{};
    }
  };


  inline match_t min_edit_distance(str_view RX, const barcode_t& bc) {
    if (RX.empty()) return { .bc = bc };  // unclear
    int min_ed = 1'000'000;
    int idx = 0;
    ReadType rt = ReadType::unclear;
    for (int i = 0; i < static_cast<int>(bc.size()); ++i) {
      auto ed = edit_distance(RX, bc[i].code(), min_ed);
      if (min_ed > ed) {
        idx = i;
        if (0 == (min_ed = ed)) [[unlikely]] {
          rt = ReadType::correct;
          break;
        }
        rt = ReadType::corrected;
      }
      else if (min_ed == ed) {
        rt = ReadType::unclear;
      }      
    }
    return { .bc = bc, .idx = idx, .ed = min_ed, .read_type = rt };
  }
  
}
