#include <iostream>
#include <iomanip>
#include <charconv>
#include <array>
#include <memory>
#include <fastq/reader.hpp>
#include <fastq/splitter.hpp>
#include <fastq/writer.hpp>


constexpr char usage_msg[] = R"(Usage: fastq_cp [OPTIONS] [FILE]
Copy ranges from fastq[.gz] files.
If FILE is not given, reads from standard input.

  -f: force overwrite of output file.
  -m <mssk>: only output unmasked lines (max. 64Bit)
    Ex: -m 0010, outputs 2nd line of every 4-line block.
  -o <FILE>: compressed output file
    If not given, writes uncompressed to standard output.
  -r <line range>: only output lines in given range.
    Ex: -r 0-10; -r 10:3
)";

// globals
std::shared_ptr<hahi::pool_t> gPool;



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


struct cin_splitter {
  std::string buf;
  
  fastq::str_view operator()() { 
    std::getline(std::cin, buf);
    return buf;
  }

  bool eof() { return std::cin.eof(); }
};


struct cout_writer {
  void puts(fastq::str_view str) {
    if (std::cin.good()) {
      std::fwrite(str.data(), 1, str.length(), stdout);
      std::fputc('\n', stdout);
    }
  }
};


template <typename Splitter, typename Writer>
void cp(Splitter&& splitter,
        Writer* writer,
        std::pair<size_t, size_t> range, 
        std::pair<uint64_t, int> mask) {
  for (size_t i = 0; !splitter.eof() && (i < range.first); ++i) {
    splitter();
  }
  auto m = mask;
  for (size_t i = range.first; !splitter.eof() && (i < range.second); ++i) {
    auto line = splitter();
    if (m.first & 1) {
      writer->puts(line);
    }
    m.first >>= 1;
    if (--m.second == 0) {
      m = mask;
    }
  }
}


int main(int argc, const char* argv[]) {
  try {
    bool force = false;
    std::pair<size_t, size_t> range{0, -1};
    std::pair<uint64_t, int> mask{-1, 64};
    std::filesystem::path file;
    std::filesystem::path output;
    int i = 1;
    while (i < argc) {
      if (0 == std::strcmp(argv[i], "-h") * std::strcmp(argv[i], "--help")) {
        throw usage_msg;
      }
      else if (0 == std::strcmp(argv[i], "-f")) {
        force = true;
      }
      else if (0 == std::strcmp(argv[i], "-r")) {
        range = parse_range((++i < argc) ? argv[i] : "");
      }
      else if (0 == std::strcmp(argv[i], "-m")) {
        mask = parse_mask((++i < argc) ? argv[i] : "");
      }
      else if (0 == std::strcmp(argv[i], "-o")) {
        if ((i + 1) < argc) {
          output = argv[++i];
        }
      }
      else if (std::filesystem::is_regular_file(argv[i])) {
        file = argv[i];
      }
      else {
        std::cerr << "invalid argument '" << argv[i] << "'\n";
        throw usage_msg;
      }
      ++i;
    }
    if (std::filesystem::exists(output) && !force) {
      throw "output file exists, consider -f";
    }
    if (!output.empty()) {
      gPool.reset( new hahi::pool_t{} );
      auto writer = std::make_unique<fastq::writer_t<>>(output, gPool);
      file.empty() ? cp(cin_splitter{}, writer.get(), range, mask)
                   : cp(fastq::line_splitter<>{file}, writer.get(), range, mask);
    }
    else {
      auto writer = std::make_unique<cout_writer>();
      file.empty() ? cp(cin_splitter{}, writer.get(), range, mask)
                   : cp(fastq::line_splitter<>{file}, writer.get(), range, mask);
    }
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
