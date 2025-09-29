#include <iostream>
#include <iomanip>
#include <charconv>
#include <array>
#include <memory>
#include <chrono>
#include <fastq/reader.hpp>
#include <fastq/splitter.hpp>
#include <fastq/writer.hpp>


constexpr char usage_msg[] = R"(Usage: fastq_paste [OPTIONS] [FILE] ...
paste line rangess from fastq[.gz] files.

  -f: force overwrite of output file.
  -m <mssk>: only output unmasked lines (max. 64Bit)
    Ex: -m 0010, outputs 2nd line of every 4-line block.
  -o <FILE>: compressed output file
    If not given, writes to standard output.
  -r <line range>: only output lines in given range.
    Ex: -r 0-10; -r 10:3
  -d: delimiter string
)";


// global thread pool
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
  if (m == 0) std::clog << "warning: empty mask\n";
  return { m, len };
}


class cout_writer {
  size_t tot_bytes_ = 0;

public:
  void puts(fastq::str_view str) {
    if (std::cin.good()) {
      std::fwrite(str.data(), 1, str.length(), stdout);
      std::fputc('\n', stdout);
      tot_bytes_ += str.length() + 1;
    }
  }
  void put(fastq::str_view str) {
    if (std::cin.good()) {
      std::fwrite(str.data(), 1, str.length(), stdout);
      tot_bytes_ += str.length();
    }
  }
  auto tot_bytes() const { return tot_bytes_; }
};


// some stats for verbose output
size_t lines_in = 0;
size_t bytes_in = 0;
size_t lines_out = 0;
size_t bytes_out = 0;


void paste(std::vector<fastq::line_splitter<>>& splitter,
           auto* writer,
           std::pair<size_t, size_t> range, 
           std::pair<uint64_t, int> mask,
           fastq::str_view delim) {
  std::vector<fastq::str_view> lines{};
  auto any_eof = [&]() { bool eof = false; for (const auto& s : splitter) { eof |= s.eof(); } return eof; };
  auto read_all = [&]() { lines.clear(); for (auto& s : splitter) { ++lines_in; lines.emplace_back(s()); } };
  for (size_t i = 0; !any_eof() && (i < range.first); ++i) {
    read_all();
  }
  auto m = mask;
  for (size_t i = range.first; !any_eof() && (i < range.second); ++i) {
    read_all();
    if (m.first & 1) {
      ++lines_out;
      writer->put(lines[0]);
      for (size_t i = 1; i < lines.size(); ++i) {       
        writer->put(delim);
        writer->put(lines[i]);
      }
      writer->put("\n");
    }
    m.first >>= 1;
    if (--m.second == 0) {
      m = mask;
    }
  }
  for (const auto& s : splitter) { bytes_in += s.tot_bytes(); };
  bytes_out += writer->tot_bytes();
}


int main(int argc, const char* argv[]) {
  try {
    // CLI arguments
    bool force = false;
    bool verbose = false;
    std::pair<size_t, size_t> range{0, -1};
    std::pair<uint64_t, int> mask{-1, 64};
    std::string delim{};
    std::vector<fastq::line_splitter<>> splitter;
    std::filesystem::path output;
    int i = 1;
    while (i < argc) {
      if (0 == std::strcmp(argv[i], "-h") * std::strcmp(argv[i], "--help")) {
        throw usage_msg;
      }
      else if (0 == std::strcmp(argv[i], "-f")) {
        force = true;
      }
      else if (0 == std::strcmp(argv[i], "-v")) {
        verbose = true;
      }
      else if (0 == std::strcmp(argv[i], "-d")) {
        if ((i + 1) < argc) {
          delim = argv[++i];
        }
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
        splitter.emplace_back(argv[i]);
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
    auto t0 = std::chrono::high_resolution_clock::now();
    if (!splitter.empty()) {
      if (!output.empty()) {
        gPool.reset( new hahi::pool_t{} );
        auto writer = std::make_unique<fastq::writer_t<>>(output, gPool);
        paste(splitter, writer.get(), range, mask, delim);
      }
      else {
        auto writer = std::make_unique<cout_writer>();
        paste(splitter, writer.get(), range, mask, delim);
      }
    }
    if (verbose) {
      auto time = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - t0);
      std::clog << "\nlines read:    " << lines_in << '\n';
      std::clog << "bytes read:    " << bytes_in << '\n';
      std::clog << "lines written: " << lines_out << '\n';
      std::clog << "bytes written: " << bytes_out << '\n';
      std::clog << "io bandwith:   " << ((bytes_in + bytes_out) / (1000*1000)) / time.count() << " MiB/s\n";
      std::clog << "elapsed time:  " << time << '\n';
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
