#include <iostream>
#include <array>
#include <fastq/reader.hpp>
#include <fastq/splitter.hpp>
#include "bench.hpp"


using splitter_t = fastq::line_splitter<>;
auto splitter = std::array<splitter_t, 4>{
  splitter_t{"../data/_gen_R1_001.fastq.gz"},
  splitter_t{"../data/_gen_R2_001.fastq.gz"},
  splitter_t{"../data/_gen_I1_001.fastq.gz"},
  splitter_t{"../data/_gen_I2_001.fastq.gz"},
};


int main() {
  try {
    const int rep = 1;
    size_t bytes_read = 0;
    size_t tot_bytes = 0;
    size_t items = 0;
    auto time = bench("", [&]() {
      while (!splitter[0].eof()) {
        for (auto& s : splitter) {
          auto view = s();
          ++items;
        }
      }
      for (const auto& s : splitter) {
        bytes_read += std::filesystem::file_size(s.reader().path());
        tot_bytes +=s.reader().tot_bytes();
      }
    }, rep);
    std::cout << "\nbytes on disk: " << bytes_read / (1000 * 1000) << " MB";
    std::cout << "\ninflated into " << double(tot_bytes) / (1000 * 1000) << " MB as " << items << " items";
    std::cout << "\ncompression rate: " << double(tot_bytes) / bytes_read << ", ";
    std::cout << double(tot_bytes / (1000 * 1000)) / time << " MB/s\n";
    return 0;
  }
  catch (const std::exception& err) {
    std::cerr << err.what() << std::endl;
  }
  return -1;
}
