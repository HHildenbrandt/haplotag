/* fastq/reader.hpp
 *
 * Copyright (c) 2025 Hanno Hildenbrandt <h.hildenbrandt@rug.nl>
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
  inline int edit_distance(str_view av, str_view bv) {
    const char* a = av.data();
    const char* b = bv.data();
    auto m = av.length();
    auto n = bv.length();
    // make outer loop shortest
    if (m > n) { std::swap(m, n); std::swap(a, b); }
    // remove matching prefixes
    while (m && (*a == *b)) { ++a; ++b; --m; --n; }
    // remove matching suffixes
    while (m && (a[m-1] == b[n-1])) { --m; --n; }
    // two lines of the distance matrix:
    alignas(64) int D0[m + 1];
    alignas(64) int D1[m + 1];
    std::iota(D0, D0 + m + 1, 0);
    auto* d0 = D0;
    auto* d1 = D1;
    for (auto i = 0; i < n; ++i) {
      const auto bi = b[i];
      d1[0] = i + 1;
      for (auto j = 0; j < m; ++j) {
        d1[j + 1] = std::min(d0[j] + (a[j] != bi), std::min(d0[j + 1], d1[j]) + 1);
      }
      std::swap(d0, d1);
    }
    return d0[m];
  }

  
  // bounded Levensthein edit distance, early exit if ed >= bound
  inline int edit_distance(str_view av, str_view bv, int bound) {
    const char* a = av.data();
    const char* b = bv.data();
    auto m = av.length();
    auto n = bv.length();
    // make outer loop shortest
    if (m > n) { std::swap(m, n); std::swap(a, b); }
    // remove matching prefixes
    while (m && (*a == *b)) { ++a; ++b; --m; --n; }
    // remove matching suffixes
    while (m && (a[m-1] == b[n-1])) { --m; --n; }
    // two lines of the distance matrix:
    alignas(64) int D0[m + 1];
    alignas(64) int D1[m + 1];
    std::iota(D0, D0 + m + 1, 0);
    auto* d0 = D0;
    auto* d1 = D1;
    for (auto i = 0; i < n; ++i) {
      const auto bi = b[i];
      d1[0] = i + 1;
      auto dmin = d0[0];
      for (auto j = 0; j < m; ++j) {
        d1[j + 1] = std::min(d0[j] + (a[j] != bi), std::min(d0[j + 1], d1[j]) + 1);
        dmin = std::min(dmin, d1[j + 1]);
      }
      if (dmin > bound) {
        return bound;   // can only get worse from here
      }
      std::swap(d0, d1);
    }
    return d0[m];
  }


  enum ReadType{
    invalid,        // code length violation
    unclear,        // multiple occurrences of same min ed
    correct,        // exact match, ed == 0
    corrected,      // unique min ed
    max_read_type
  };


  struct match_t {
    int idx = 0;          // index into bc, 0 if read_type == unclear
    int ed = -1;          // edit distance
    ReadType rt = ReadType::invalid;
  };


  inline match_t min_edit_distance(str_view RX, size_t code_length, const barcode_t& bc) {
    if (RX.length() < code_length) return {};  // invalid
    int min_ed = 1'000'000;
    int idx = 0;
    ReadType rt = ReadType::unclear;
    for (size_t i = 1; i < bc.size(); ++i) {
      auto ed = edit_distance(RX, bc[i].code, min_ed + 1);
      if (min_ed > ed) {
        idx = i;
        if (0 == (min_ed = ed)) [[unlikely]] {
          rt = ReadType::correct;
          break;  // assuming unique barcodes
        }
        rt = ReadType::corrected;
      }
      else if (min_ed == ed) {
        rt = ReadType::unclear;
      }      
    }
    return { .idx = (rt == ReadType::unclear) ? 0 : idx, .ed = min_ed, .rt = rt };
  }
  

  // min_edit_distance, assuming code_length = bc.max_code_length()
  inline match_t min_edit_distance(str_view RX, const barcode_t& bc) {
    return min_edit_distance(RX, bc.max_code_length(), bc);
  }
  
}
