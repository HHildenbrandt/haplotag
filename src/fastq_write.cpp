#include <iostream>
#include <random>
#include <array>
#include <fastq/writer.hpp>
#include "bench.hpp"


using namespace std::literals::chrono_literals;


std::string line = "jgfaj;ltorglkn,.zdfgmdlkgb.c .zxmglfdk;bmbz/.gkrdxjg;lfkb";
constexpr int npermute = 1;
auto reng = std::mt19937{};
auto ud = std::uniform_int_distribution<size_t>(4, line.size() - 1);


void permute() {
  for (auto p = 0; p < npermute; ++p) {
    std::swap(line[ud(reng)], line[ud(reng)]);
  };
}


int main(int argc, char** argv) {
  size_t MB = 500 * 1000 * 1000;
  size_t tot_bytes = 0;
  unsigned nt = 4;
  if (argc > 1) {
    nt = std::atoi(argv[1]);
  }
  std::shared_ptr<hahi::pool_t> gPool{ new hahi::pool_t{nt} };
  size_t lines = MB / line.size();
  lines = 4 * ((lines + 3) / 4);    // align to fastq sequence
  auto writers = std::array<fastq::writer_t<>, 4>{
    fastq::writer_t(std::string("../data/_dummy_writer") + "0.gz", gPool),
    fastq::writer_t(std::string("../data/_dummy_writer") + "1.gz", gPool),
    fastq::writer_t(std::string("../data/_dummy_writer") + "2.gz", gPool),
    fastq::writer_t(std::string("../data/_dummy_writer") + "3.gz", gPool),
  };
  std::cout << "Using " << gPool->num_threads() << " compression thread(s)" << std::endl;
  auto time = bench("", [&]() {
    for (size_t i = 0; i < lines; ++i) {
      permute();
      for (auto& writer_t : writers) {
        writer_t.puts(line);
      }
    }
    // be fair
    for (auto& writer_t : writers) {
      writer_t.close(true);
      tot_bytes += writer_t.tot_bytes();
    }
  }, 1);
  std::cout << "  ~ " << (tot_bytes / MB) / time << " MB/s" << std::endl;
  return 0;
}
