#include <cstring>
#include <iostream>
#include <fastq/reader.hpp>
#include <fastq/splitter.hpp>


using splitter_t = fastq::line_splitter<>;


int main() {
  auto s = splitter_t("../data/R1_001.fastq.gz");
  size_t items = 0;
  while (!s.eof()) {
    auto x = s(4);
    ++items;
  }
  std::cout << items << " items read, " << s.reader().tot_bytes() / (1000 * 1000) << " MB" << std::endl;
  return 0;
}