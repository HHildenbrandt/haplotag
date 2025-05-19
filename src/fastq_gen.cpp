#include <iostream>
#include <random>
#include <list>
#include <algorithm>
#include <fastq/splitter.hpp>
#include "fastq/reader.hpp"
#include "fastq/writer.hpp"


auto reng = std::mt19937{};

void permute(std::string_view str) {
  if (str.length() <4) return;
  auto dist = std::uniform_int_distribution<size_t>(1, str.length() - 2);
  auto ptr = const_cast<char*>(str.data()); // ugly as hell
  std::swap(*(ptr + dist(reng)), *(ptr + dist(reng)));
}


int main(int argc, const char* argv[]) {
  size_t N = 1;    // duplicates
  if (argc > 1) N = std::stoull(argv[1]);

  // shared thread-pool
  std::shared_ptr<hahi::pool_t> shPool{ new hahi::pool_t(-1) };
  
  using splitter_t = fastq::seq_field_splitter<>;
  auto splitter = std::list<splitter_t>{};
  auto writer = std::list<fastq::writer_t<>>{};
  auto data_dir = std::filesystem::path("../data/");
  for (auto const& dir_entry : std::filesystem::directory_iterator{data_dir}) {
    if (dir_entry.is_regular_file()) {
      auto p = dir_entry.path();
      auto file = p.filename();
      if ((file.extension() == ".gz") && !file.string().starts_with("_")) {
        std::cout << "found '" << p.string() << '\'' << std::endl;
        splitter.emplace_back(p);
        writer.emplace_back(p.replace_filename(std::string("_gen_") + file.string()), shPool);
      }
    }
  }
  std::cout << "generating " << writer.size() << " '_gen_*.gz' files..." << std::endl;
  while (!std::all_of(splitter.begin(), splitter.end(), [](auto& s) { return s.eof(); })) {
    auto wit = writer.begin();
    for (auto sit = splitter.begin(); sit != splitter.end(); ++sit, ++wit) {
      if (!sit->eof()) {
        auto seq = (*sit)();
        for (size_t n = 0; n < N; ++n) {
          for (auto& field : seq) {
            if (N > 1) permute(field);
            wit->puts(field);
          }
        }
      }
    }
  }
  size_t tot_bytes_inf = 0;
  size_t tot_bytes_def = 0;
  auto sit = splitter.begin();
  auto wit = writer.begin();
  for (; sit != splitter.end(); ++sit, ++wit) {
    wit->close(true);
    tot_bytes_inf += sit->reader().tot_bytes();
    tot_bytes_def += wit->tot_bytes();
  }
  std::cout << double(tot_bytes_inf) / (1000 * 1000) << " MB inflated" << std::endl;
  std::cout << double(tot_bytes_def) / (1000 * 1000 * 1000) << " GB deflated" << std::endl;
  return 0;
}
