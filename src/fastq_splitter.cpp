#include <cstring>
#include <iostream>
#include <fastq/reader.hpp>
#include <fastq/splitter.hpp>


using reader_t = fastq::reader_t;
using splitter_t = fastq::seq_splitter<reader_t>;


int main() {
  auto s = splitter_t("../data/_gen_I2_001.fastq.gz");
  size_t items = 0;
  while (!s.eof()) {
    auto x = s();
    ++items;
  }
  auto tot_bytes = s.reader().tot_bytes();
  std::cout << items << " items read, " << tot_bytes / (1000 * 1000) << " MB" << std::endl;
  return 0;
}