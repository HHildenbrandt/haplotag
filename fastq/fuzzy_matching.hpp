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
  inline size_t edit_distance(str_view av, str_view bv) {
    const char* a = av.data();
    const char* b = bv.data();
    auto m = av.length();
    auto n = bv.length();
    if (m > n) { std::swap(m, n); std::swap(a, b); }
    // remove matching prefixes and suffixes
    while (m && (*a == *b)) { ++a; ++b; --m; --n; }
    while (m && (a[m-1] == b[n-1])) { --m; --n; }
    size_t D[m + 1];  // single row of the distance matrix
    std::iota(D, D + m + 1, 0);
    for (auto i = 1; i <= n; ++i) {
      const auto bi = b[i - 1];
      auto tmp = std::exchange(D[0], i);  // tmp <- L(i-1,0), L(i,0) <- i
      for (auto j = 1; j <= m; ++j) {
        if (a[j - 1] != bi) {
          tmp = std::min(D[j], std::min(D[j - 1], tmp)) + 1;
        }
        std::swap(tmp, D[j]);   // d[j] <- L(i,j)
      }
    }
    return D[m];
  }
  

  // bounded Levensthein edit distance, early exit if ed >= bound
  inline size_t edit_distance(str_view av, str_view bv, size_t bound) {
    const char* a = av.data();
    const char* b = bv.data();
    auto m = av.length();
    auto n = bv.length();
    if (m > n) { std::swap(m, n); std::swap(a, b); }
    // remove matching prefixes and suffixes
    while (m && (*a == *b)) { ++a; ++b; --m; --n; }
    while (m && (a[m-1] == b[n-1])) { --m; --n; }
    size_t D[m + 1];  // single row of the distance matrix
    std::iota(D, D + m + 1, 0);
    for (auto i = 1; i <= n; ++i) {
      const auto bi = b[i - 1];
      auto dmin = D[0];
      auto tmp = std::exchange(D[0], i);  // tmp <- L(i-1,0), L(i,0) <- i
      for (auto j = 1; j <= m; ++j) {
        if (a[j - 1] != bi) {
          tmp = std::min(D[j], std::min(D[j - 1], tmp)) + 1;
        }
        dmin = std::min(dmin, tmp);
        std::swap(tmp, D[j]);   // d[j] <- L(i,j)
      }
      if (dmin >= bound) {
        return bound;    // can only get worse...
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
    int idx = 0;          // index into bc, 0 if read_type == unclear
    int ed = -1;          // edit distance
    ReadType read_type = ReadType::unclear;
    
    str_view tag() const { return bc[idx].tag; }
    str_view code() const { return bc[idx].code; }
  };


  inline match_t min_edit_distance(str_view RX, const barcode_t& bc) {
    if (RX.empty()) return { .bc = bc };  // unclear
    int min_ed = 1'000'000;
    int idx = 0;
    ReadType rt = ReadType::unclear;
    for (int i = 1; 1 < static_cast<int>(bc.size()); ++i) {
      auto ed = edit_distance(RX, bc[i].code, min_ed + 1);
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
    return { .bc = bc, .idx = (rt == ReadType::unclear) ? 0 : idx, .ed = min_ed, .read_type = rt };
  }
  
}
