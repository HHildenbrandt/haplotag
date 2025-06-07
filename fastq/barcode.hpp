/* fastq/barcode.hpp
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
#include <stdexcept>
#include <fstream>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>
#include "fastq.hpp"
#include "splitter.hpp"


namespace fastq {


  // returns empty string_view if [pos, pos+count) is out of bounds
  inline str_view save_substr(str_view str, 
              size_t pos,
              size_t count = str_view::npos) noexcept {
    if (pos > str.length()) [[unlikely]] return {};
    auto rlen = std::min(str.length() - pos, count);
    if ((count != str_view::npos) && (rlen < count)) return {};
    auto x = str_view{ str.data() + pos, rlen };
    return x;
  }
  

  // barcode file reader
  // more general: <tag, code> file reader
  class barcode_t {
  public:
    struct entry_t {
      explicit entry_t(const std::string& line)
      : code_(line.substr(line.find_last_of("\t ") + 1)), // tab or space
        tag_(line.substr(0, line.find_first_of("\t ")))   // tab or space
      {
        if (tag_.empty() || code_.empty()) {
          throw std::runtime_error("corrupted barcode line");
        }
      }

      char code_letter() const noexcept { return tag_[0]; }
      const std::string& code() const noexcept { return code_; }
      const std::string& tag() const noexcept { return tag_; }

    private:
      friend class barcode_t;
      const std::string code_;
      const std::string tag_;
    };

    // barcode_t(barcode_t&) = delete;
    // barcode_t& operator=(barcode_t&) = delete;
    barcode_t() {}

    explicit barcode_t(const std::filesystem::path& path, std::string unclear_tag = {}) : path_(path) {
      try {
        auto is = std::ifstream(path);
        auto lines = std::vector<std::string>{};
        std::string line{};
        while (is) {
          std::getline(is, line);
          if (!line.empty()) lines.push_back(std::move(line));  // last line allowed to be '\n'
        }
        if (lines.empty()) throw std::runtime_error("doesn't exist or corrupted");
        const auto probe = entry_t(lines[0]);
        if (unclear_tag.empty()) {
          unclear_tag = std::string{probe.code_letter()} + std::string(probe.tag().length() - 1, '0');
        }
        unclear_tag_ = unclear_tag;
        for (const auto& line : lines) {
          auto& bc = bc_.emplace_back(line);
          if (bc.tag() == unclear_tag) {
            throw std::runtime_error("duplicated dummy-tag.");
          }
          min_code_length_ = std::min(min_code_length_, bc.code().length());
          max_code_length_ = std::max(max_code_length_, bc.code().length());
        }
      }
      catch (const std::exception& err) {
        throw std::runtime_error(path.string() + ": " + err.what());
      }
    }

    bool empty() const noexcept { return bc_.empty();}
    size_t size() const noexcept { return bc_.size();}
    size_t min_code_length() const noexcept { return min_code_length_; }
    size_t max_code_length() const noexcept { return max_code_length_; }
    
    const entry_t& operator[](size_t idx) const { return bc_[idx];}
    auto begin() const { return bc_.begin(); }
    auto cbegin() const { return bc_.cbegin(); }
    auto end() const { return bc_.end(); }
    auto cend() const { return bc_.cend(); }
    auto data() const { return bc_.data(); }
    
    const std::string& unclear_tag() const noexcept { return unclear_tag_; }
    const std::filesystem::path& path() const noexcept { return path_; }

  private:
    std::vector<entry_t> bc_;
    std::string unclear_tag_;
    size_t min_code_length_ = 1'000'000;
    size_t max_code_length_ = 0;
    std::filesystem::path path_;
  };


}

