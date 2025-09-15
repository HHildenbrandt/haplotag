/*
* Copyright (c) 2019-2025 Andreea Dreau https://github.com/adreau
*                         Frank Chan https://github.com/evolgenomics
*                         Hanno Hildenbrandt <h.hildenbrandt@rug.nl>
*                         Thijs Janzen
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


#include <iostream>
#include <filesystem>
#include <numeric>
#include <deque>
#include <map>
#include <future>
#include <thread>
#include <type_traits>
#include "device/pool.hpp"
#include "fastq/barcode.hpp"
#include <fastq/reader.hpp>
#include "fastq/splitter.hpp"
#include "fastq/writer.hpp"
#include "fastq/fuzzy_matching.hpp"

// uncomment to reproduce match-collapsing in demult_fastq
#define DEMULT_DASTQ_BEHAVIOUR


using namespace fastq;
namespace fs = std::filesystem;
using namespace std::chrono_literals;

// global thread-pool
std::shared_ptr<hahi::pool_t> gPool{ new hahi::pool_t{unsigned(-1)} };

const auto bc_A = barcode_t("../data/BC_A.txt");
const auto bc_B = barcode_t("../data/BC_B.txt");
const auto bc_C = barcode_t("../data/BC_C.txt");
const auto bc_D = barcode_t("../data/BC_D.txt");

using Rsplitter = line_splitter<>;
using Isplitter = seq_field_splitter<0b1010>;   // [RX, QX]

auto R1_in = Rsplitter("../data/_gen_R1_001.fastq.gz");
auto R2_in = Rsplitter("../data/_gen_R2_001.fastq.gz");
auto I1_in = Isplitter("../data/_gen_I1_001.fastq.gz");
auto I2_in = Isplitter("../data/_gen_I2_001.fastq.gz");

enum MaskedNames {
  RName = 0,
  IRX = 0,
  IQX = 1
};

auto R1_out = writer_t<>("../data/_fuzzy_R1_001.fastq.gz", gPool);
auto R2_out = writer_t<>("../data/_fuzzy_R2_001.fastq.gz", gPool);

// log files
auto clear_log = fs::path{"../data/_fuzzy_clearBC.log"};
auto unclear_log = fs::path{"../data/_fuzzy_unclearBC.log"};


// zero-initialized std::array
struct match_counts : std::array<size_t, max_read_type> {
  match_counts() { this->fill(0); }
};
using read_type_map_t = std::map<std::string, match_counts>;  // <code, counters>


struct blk_fuzzy_out_t {
  read_type_map_t rm;
  std::vector<std::string> names;
};


// runs on multiple threads
blk_fuzzy_out_t blk_fuzzy(Isplitter::blk_type&& I1, Isplitter::blk_type&& I2) {
  assert(I1.size() == I2.size());
  std::string code{};
  read_type_map_t read_map{};
  std::vector<std::string> names(I1.size());
  for (size_t i = 0; i < I1.size(); ++i) {
    auto I1_RX = I1[i][IRX];
    auto I2_RX = I2[i][IRX];
    match_t matches[] = {
      min_edit_distance(max_substr(I1_RX, 7, 6), bc_A),
      min_edit_distance(max_substr(I1_RX, 0, 6), bc_C),
      min_edit_distance(max_substr(I2_RX, 7, 6), bc_B),
      min_edit_distance(max_substr(I2_RX, 0, 6), bc_D)
    };
#if defined(DEMULT_DASTQ_BEHAVIOUR)
    if ((matches[0].read_type == unclear) || (matches[1].read_type == unclear)) {
      matches[0].read_type = matches[1].read_type = unclear;
    }
#endif
    code.clear();
    // accumulate read_types counts
    auto acc = match_counts{};
    for (size_t j = 0; j < std::size(matches); ++j) {
      code.append(matches[j].tag());
      for (int s = 0; s < max_read_type; ++s) {
        acc[s] += (matches[j].read_type == s);
      }
    }
    // reduction
    ReadType match = (4 == acc[correct]) 
                   ? correct
                   : ((0 < acc[unclear]) ? unclear : corrected);
    ++read_map[code][match];   // map[key] inserts <key, {}> or retrieves <key, val>
    // assemble new 'name' field suffix
    names[i].append("BX:Z:").append(code);
    names[i].append("\tRX:Z:").append(I1_RX).append("+").append(I2_RX);
    names[i].append("\tQX:Z:").append(I1[i][IQX]).append("+").append(I2[i][IQX]);
  }
  return { std::move(read_map), std::move(names) };
}


// returns number of processes items
// runs on main thread
size_t merge_blk_fuzzy(read_type_map_t& dst, blk_fuzzy_out_t&& bf) {
  for (const auto& [code, c] : bf.rm) {
    for (auto i = 0; i < max_read_type; ++i) {
      dst[code][i] += c[i];
    }
  }
  // generate output, try to avoid (re)allocations
  for (const auto& name : bf.names) {
    auto R1 = R1_in();
    R2_in();  // skip
    R1 = R1.substr(RName, R1.find_first_of("\t ") + 1);   // tab or space
    R1_out.put(R1); R1_out.puts(name);
    R2_out.put(R1); R2_out.puts(name);  // copy over
    for (int i = 1; i < 4; ++i) {
      R1_out.puts(R1_in());
      R2_out.puts(R2_in());
    }
  }
  return bf.names.size();
}


int main(int argc, const char * argv[]) {
  static constexpr size_t max_block_size = 100'000;
  size_t block_size = 10'000;
  size_t csets = 0;
  try {
    std::deque<std::future<blk_fuzzy_out_t>> futures;
    auto read_type_map = read_type_map_t{};
    std::cout << "starting" << std::endl;
    while (!I1_in.eof()) {
      assert(!I2_in.eof());
      futures.emplace_back(
        gPool->async(blk_fuzzy, I1_in(block_size), I2_in(block_size))
      );
      while (!futures.empty() && (std::future_status::ready == futures.front().wait_for(0s))) {
        csets += merge_blk_fuzzy(read_type_map, futures.front().get());
        futures.pop_front();
        std::cout << "*  " << csets << " sets processed +\r" << std::flush;
        block_size = std::min((4 * block_size) / 3, max_block_size);
      }
    }
    // process results from remaining blocks
    while (!futures.empty()) {
      csets += merge_blk_fuzzy(read_type_map, futures.front().get());
      futures.pop_front();
    }
    std::cout << "*  " << csets << " sets processed -\r" << std::flush;

    // split log map into two log files
    auto os_clear = std::ofstream(clear_log);
    auto os_unclear = std::ofstream(unclear_log);
    os_clear << "Barcode \t Correct reads \t Corrected reads\n";
    os_unclear << "Barcode \t Reads\n";
    for (const auto& [code, cm] : read_type_map) {
      if (cm[correct] || cm[corrected]) {
        os_clear << code << '\t' << cm[correct] << '\t' << cm[corrected] << '\n';
      }
      else {
        os_unclear << code << '\t' << cm[unclear] << '\n';
      }
    }
    R1_out.close(true);   // early close & flush to make failed flag and tot_bytes visible
    R2_out.close(true);   // early close & flush to make failed flag and tot_bytes visible
    if (R1_in.failed() || R2_in.failed() || I1_in.failed() || I2_in.failed() || R1_out.failed() || R2_out.failed()) {
      throw std::runtime_error("something went horrible wrong with the gz files");
    }
    const size_t tb_com = R1_out.tot_bytes() + R2_out.tot_bytes();
    const size_t tb_dec = R1_in.reader().tot_bytes() + R2_in.reader().tot_bytes() + I1_in.reader().tot_bytes() + I2_in.reader().tot_bytes();
    std::cout << "\n*  " << tb_dec / (1000 * 1000) << " MB decompressed" << std::endl;
    std::cout << "*  " << tb_com / (1000 * 1000) << " MB compressed" << std::endl;
    return 0;
  }
  catch (const std::exception& err) {
    std::cerr << err.what() << std::endl;
  }
  catch (...) {
    std::cerr << "Unknown exception" << std::endl;
  }
  return -1;
}
