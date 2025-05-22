#include <cstring>
#include <iostream>
#include <fastq/reader.hpp>
#include <fastq/splitter.hpp>


using splitter_t = fastq::seq_splitter<>;


int main() {
  auto s = splitter_t("../data/_gen_test_H4_R2_001.fastq.gz");
  size_t items = 0;
  while (!s.eof()) {
    auto x = s();
    if (s.eof()) {
      int dummy = 0;
    }
    ++items;
  }
  auto tot_bytes = s.reader().tot_bytes();
  std::cout << items << " items read, " << tot_bytes / (1000 * 1000) << " MB" << std::endl;
  return 0;
}