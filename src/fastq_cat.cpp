#include <iostream>
#include <charconv>
#include <array>
#include <fastq/reader.hpp>
#include <fastq/splitter.hpp>


constexpr char usage_msg[] = R"(Usage: fastq_cat [OPTION] [FILE]

  -m <mssk>: only output unmasked lines
    Ex: -m 1001
  -n: show line numbers
  -r <line range>: only output lines in given range
    Ex: -r 0-10; -r 10:3
  -t: output tail 
)";


std::pair<size_t, size_t> parse_range(std::string_view str) {
  size_t n0 = 0;
  size_t n1 = size_t(-1);
  if (str.empty()) return { n0, n1 };
  auto [p0, ec0] = std::from_chars(str.begin(), str.end(), n0);
  if (ec0 != std::errc{}) throw "can't parse range";
  if (p0 == str.end()) return { n0, n1 };
  char delim = *p0++;
  if (!((delim == '-') || (delim == ':'))) throw "can't parse range";  
  auto [p1, ec1] = std::from_chars(p0, str.end(), n1);
  if (ec1 != std::errc{}) throw "can't parse range\n";
  return { n0, (delim == ':') ? n0 + n1 : n1 };
}


std::pair<uint64_t, int> parse_mask(std::string_view str) {
  uint64_t m = 0;
  int len = static_cast<int>(str.length());
  if (len > 64) throw "mask exceeds 64 bit";
  for (auto chr : str) {
    m <<= 1;
    if (chr == '1') m |= 1;
    else if (chr != '0') throw "can't parse mask";
  }
  if (m == 0) throw "empty mask";
  return { m, len };
}


void cat(std::pair<size_t, size_t> range, std::pair<uint64_t, int> mask, const std::filesystem::path& file) {
  auto splitter = fastq::line_splitter<>(file);
  for (size_t i = 0; !splitter.eof() && (i < range.first); ++i) {
    splitter();
  }
  auto m = mask;
  for (size_t i = range.first; !splitter.eof() && (i < range.second); ++i) {
    auto line = splitter();
    if (m.first & 1) {
      std::cout << line << '\n';
    }
    m.first >>= 1;
    if (--m.second == 0) {
      m = mask;
    }
  }
}


int main(int argc, const char* argv[]) {
  auto x = parse_mask("1011");
  try {
    bool tail = false;
    std::pair<size_t, size_t> range{0, -1};
    std::pair<uint64_t, int> mask{-1, 64};
    std::filesystem::path file;
    int i = 1;
    while (i < argc) {
      if (0 == std::strcmp(argv[i], "-h") * std::strcmp(argv[i], "--help")) {
        throw usage_msg;
      }
      if (0 == std::strcmp(argv[i], "-t")) {
        tail = true;
      }
      else if (0 == std::strcmp(argv[i], "-r")) {
        range = parse_range((++i < argc) ? argv[i] : "");
      }
      else if (0 == std::strcmp(argv[i], "-m")) {
        mask = parse_mask((++i < argc) ? argv[i] : "");
      }
      else if (std::filesystem::is_regular_file(argv[i])) {
        file = argv[i];
      }
      else throw "unknown argument";
      ++i;
    }
    if (file.empty()) {
      throw "file not given or inaccessible";
    }
    cat(range, mask, file);
    return 0;
  }
  catch (const char* err) {
    std::cerr << err << '\n';
  }
  catch (const std::exception& err) {
    std::cerr << err.what() << '\n';
  }
  return -1;
}
