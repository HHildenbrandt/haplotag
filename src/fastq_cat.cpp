#include <iostream>
#include <charconv>
#include <array>
#include <fastq/reader.hpp>
#include <fastq/splitter.hpp>


using splitter_t = fastq::line_splitter<>;


bool is_gz_file(const std::filesystem::path& path) {
  return (".gz" == path.extension()) && std::filesystem::is_regular_file(path);
}


std::pair<size_t, size_t> range(std::string_view range_str) {
  size_t from = 0;
  size_t to = 0;
  auto [ptr, ec] = std::from_chars(range_str.begin(), range_str.end(), from);
  if (ec != std::errc{}) throw "can't parse range argument";
  if (ptr != range_str.end()) {
    auto [_, ec] = std::from_chars(++ptr, range_str.end(), to);
    if (ec != std::errc{}) throw "can't parse range argument";
  }
  else {
    std::swap(from, to);
  }
  return {from, to};
}

int main(int argc, const char* argv[]) {
  try {
    bool tail = false;
    bool seq = false;
    std::filesystem::path input;
    std::string range_str{};
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
    if (input.empty()) {
      throw "input not given or inaccessible";
    }
    std::cout << tail << '\n';
    std::cout << seq << '\n';
    std::cout << input << '\n';
    std::cout << range_str << '\n';
    auto x = range(range_str);
    return 0;
  }
  catch (const char* err) {
    std::cerr << err << std::endl;
  }
  catch (const std::exception& err) {
    std::cerr << err.what() << std::endl;
  }
  return -1;
}
