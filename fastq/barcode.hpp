/* fastq/barcode.hpp
 *
 * Copyright (c) 2025 Hanno Hildenbrandt <h.hildenbrandt@rug.nl>
 */
 
#pragma once

#include <cassert>
#include <stdexcept>
#include <fstream>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>
#include <algorithm>
#include "fastq.hpp"
#include "splitter.hpp"


namespace fastq {


  // barcode file reader
  // more general: <tag, code> file reader
  class barcode_t {
  public:
    struct entry_t {
      std::string tag;
      std::string code;
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
          // last line allowed to be '\n' 
          if (!line.empty()){ 
            lines.push_back(std::move(line));
          }
        }
        if (lines.empty()) throw std::runtime_error("doesn't exist or is corrupted");
        bc_.emplace_back();   // placeholder unclear
        for (const auto& line : lines) {
          const auto p0 = line.find_first_of(" \t");  // space or tab
          if (p0 == line.npos) throw std::runtime_error("corrupted line");
          auto& bc = bc_.emplace_back(entry_t{.tag = line.substr(0, p0), .code = line.substr(p0 + 1)});
          min_code_length_ = std::min(min_code_length_, bc.code.length());
          max_code_length_ = std::max(max_code_length_, bc.code.length());
        }
        if (unclear_tag.empty()) {
          // create heuristic 'code_letter' 0+ tag
          unclear_tag = std::string(bc_[1].tag.size(), '0');
          unclear_tag[0] = bc_[1].tag[0];
        }
        bc_[0].tag = unclear_tag;
      }
      catch (const std::exception& err) {
        throw std::runtime_error(path.string() + ": " + err.what());
      }
    }

    // replace first tag-character with code_letter
    void reset_code_letter(char code_letter) {
      for (auto& e : bc_) {
        e.tag[0] = code_letter;
      }
    }

    // sort true entries by tag
    void sort_by_tags() {
      if (bc_.size() > 1) {
        std::sort(bc_.begin() + 1, bc_.end(), [](const auto& a, const auto& b) {
          return a.tag < b.tag;
        });
      }
    }

    bool empty() const noexcept { return bc_.empty();}
    size_t size() const noexcept { return bc_.size();}
    size_t min_code_length() const noexcept { return min_code_length_; }
    size_t max_code_length() const noexcept { return max_code_length_; }
    
    // idx == 0 returns unclear_tag
    const entry_t& operator[](size_t idx) const { return bc_[idx];}
    auto begin() const { return bc_.begin(); }
    auto cbegin() const { return bc_.cbegin(); }
    auto end() const { return bc_.end(); }
    auto cend() const { return bc_.cend(); }
    auto data() const { return bc_.data(); }
    
    //const std::string& unclear_tag() const noexcept { return unclear_tag_; }
    const std::filesystem::path& path() const noexcept { return path_; }

  private:
    std::vector<entry_t> bc_;
    size_t min_code_length_ = 1'000'000;
    size_t max_code_length_ = 0;
    std::filesystem::path path_;
  };


}

