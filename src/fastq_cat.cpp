#include <iostream>
#include <array>
#include <fastq/reader.hpp>
#include <fastq/splitter.hpp>


using splitter_t = fastq::line_splitter<>;


bool is_gz_file(const std::filesystem::path& path) {
  return (".gz" == path.extension()) && std::filesystem::is_regular_file(path);
}


int main(int argc, const char* argv[]) {
  try {
    bool tail = false;
    bool seq = false;
    std::filesystem::path input;
    std::string range_str{"0:10"};
    int i = 1;
    while (i < argc) {
      if (0 == std::strcmp(argv[i], "-t")) {
        tail = true;
        ++i;
        continue;
      }
      if (0 == std::strcmp(argv[i], "-s")) {
        seq = true;
        ++i;
        continue;
      }
      if (is_gz_file(argv[i])) {
        input = argv[i++];
        continue;
      }
      range_str = argv[i++];
    }
    std::cout << tail << '\n';
    std::cout << seq << '\n';
    std::cout << input << '\n';
    std::cout << range_str << '\n';
    return 0;
  }
  catch (const std::exception& err) {
    std::cerr << err.what() << std::endl;
  }
  return -1;
}
