#include <iostream>
#include <filesystem>
#include <algorithm>
#include <nlohmann/json.hpp>
#include <fastq/barcode.hpp>
#include <fastq/reader.hpp>
#include <fastq/splitter.hpp>
#include <fastq/writer.hpp>
#include <fastq/fuzzy_matching.hpp>
#include "device/pool.hpp"


constexpr char usage_msg[] = R"(Usage: fastq_H4 JSON_FILE [OPTIONS]...
  -h, --help: show this message
  -f, --force: force overwrite of output directory.
  --replace "json expr"
    Ex: --replace '"/range:" "0-1000"'
  --dry: dry-run
)";


namespace fs = std::filesystem;
using namespace std::chrono_literals;
using json = nlohmann::json;


// global thread-pool
std::shared_ptr<hahi::pool_t> gPool;
using Splitter = fastq::seq_field_splitter<0b1111>;


// parse "range" value
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


struct H4 {
  static constexpr size_t blk_size = 10000;   // tuneable

  explicit H4(const nlohmann::json& J);

  bool has_stagger() const noexcept { return !stagger.empty(); }
  bool has_plate() const noexcept { return !plate.empty(); }

  void dry_run();
  void run();

  std::pair<size_t, size_t> range;
  fastq::barcode_t bc_A;
  fastq::barcode_t bc_B;
  fastq::barcode_t bc_C;
  fastq::barcode_t bc_D;
  fastq::barcode_t plate;
  fastq::barcode_t stagger;

  fastq::seq_field_splitter<0b1111> R1;
  Splitter R2;
  Splitter R3;
  Splitter R4;
  Splitter I1;

  std::unique_ptr<fastq::writer_t<>> R1_out;
  std::unique_ptr<fastq::writer_t<>> R2_out;

  bool verbose = false;
  bool clipping = false;
  std::filesystem::path bc_root;
  std::filesystem::path gz_sub;
  std::filesystem::path out_dir;
};


#define optional_json(expr) try { expr; } catch (json::exception&) {}


H4::H4(const json& J) {
  auto root = J.at("root").get<fs::path>();
  range = parse_range(J.at("range").get<std::string>());
  if (range.first >= range.second) throw "invalid range";
  gPool.reset( new hahi::pool_t(J.at("num_threads").get<unsigned>()));

  // barcodes
  auto jbc = J.at("barcodes");
  auto bc_root = root;  // default: same as /root
  if (auto bc_root_str = jbc.at("root").get<std::string>(); !bc_root_str.empty()) {
    bc_root = bc_root_str;  // overwrite if not empty
  }
  auto gen_bc = [&](const char* L) {
    auto bc = fastq::barcode_t{bc_root / jbc.at(L).at("file").get<std::string>(), jbc.at(L).at("unclear_tag")};
    optional_json(bc.reset_code_letter(jbc.at(L).at("code_letter").get<std::string>()[0]));
    bool sort = false;
    optional_json(sort = jbc.at(L).at("sort_by_tag").get<bool>());
    if (sort) bc.sort_by_tags();
    return bc;
  };
  bc_A = gen_bc("A");
  bc_B = gen_bc("B");
  bc_C = gen_bc("C");
  bc_D = gen_bc("D");
  optional_json(plate = gen_bc("plate"));
  optional_json(stagger = gen_bc("stagger"));

  // reads
  auto jr = J.at("reads");
  gz_sub = root.string() + jr.at("subdir").get<std::string>();
  R1 = Splitter{gz_sub / jr.at("R1").get<std::string>()};
  R2 = Splitter{gz_sub / jr.at("R2").get<std::string>()};
  R3 = Splitter{gz_sub / jr.at("R3").get<std::string>()};
  R4 = Splitter{gz_sub / jr.at("R4").get<std::string>()};
  I1 = Splitter{gz_sub / jr.at("I1").get<std::string>()};
  
  // output
  auto jout = J.at("output");    
  out_dir = root / jout.at("subdir").get<std::string>();
  R1_out.reset(new fastq::writer_t<>{out_dir / jout.at("R1").get<std::string>(), gPool});
  R2_out.reset(new fastq::writer_t<>{out_dir / jout.at("R2").get<std::string>(), gPool});
}


void H4::run() {
  using matches_t = std::vector<fastq::match_t>;

  size_t i = 0;
  // skip to range.first
  for (; i < range.first; ++i) {
    for (auto* R : { &R1, &R2, &R3, &R4, &I1 }) {
      if (R->eof()) break;
      R->operator()();
    }
  }
  if (i != range.first) throw "range exceeds number of reads";
  // work pipeline
  auto future_queue = std::deque<std::future<matches_t>>{};
  bool any_eof = false;
  for (; !any_eof && (i < range.second); i += blk_size) {
    // collect max. blk_size sequences from all fastq.gz files
    auto blks = std::vector<Splitter::blk_type>{};
    for (auto* R : { &R1, &R2, &R3, &R4, &I1 }) {
      any_eof |= R->eof();
      blks.emplace_back(R->operator()(blk_size));
    }
    future_queue.emplace_back(gPool->async([this, blks = std::move(blks)]() {
      auto matches = matches_t{};
      matches.reserve(blk_size);
      for (auto& seq : blks[3]) {
        matches.push_back(fastq::min_edit_distance(fastq::max_substr(seq[1], 0, 10), bc_A));
      }
      return matches;
    }));
    while (!future_queue.empty() && (std::future_status::ready == future_queue.front().wait_for(0s))) {
      auto m = std::move(future_queue.front().get());
      future_queue.pop_front();
    }
  }
  while (!future_queue.empty()) {
    auto m = std::move(future_queue.front().get());
    future_queue.pop_front();
  }
  int dummy = 0;
}


void H4::dry_run() {
  using std::cout;
  cout << "range: " << range.first << '-' << range.second << '\n';
  cout << "num_threads: " << gPool->num_threads() << '\n';
  auto bc_stats = [](const char* name, const auto& bc) { 
    cout << name << "  ";
    if (bc.empty()) {
      cout << "NA\n";
      return;
    }
    cout << '"' << bc[0].tag << "\"  "
         << bc.size() -1 << "  "
         << '[' << bc.min_code_length() << ", " << bc.max_code_length() << "]  "
         << bc.path()
         << '\n';
  };
  cout << "barcodes\n";
  bc_stats("    bc_A:   ", bc_A);
  bc_stats("    bc_B:   ", bc_B);
  bc_stats("    bc_C:   ", bc_C);
  bc_stats("    bc_D:   ", bc_D);
  bc_stats("    plate:  ", plate);
  bc_stats("    stagger:", stagger);

  auto gz_stats = [](const char* name, const auto& gz) {
    cout << name << "  ";
    if (gz.failed()) cout << "NA\n";
    else cout << gz.reader().path() << '\n';
  };
  cout << "reads\n";
  gz_stats("    R1:", R1);
  gz_stats("    R2:", R2);
  gz_stats("    R3:", R3);
  gz_stats("    R4:", R4);
  gz_stats("    I1:", I1);

  cout << "matches\n";
  if (has_stagger()) {
    std::cout << "    stagger <- idx min_ed(R2[1](0:" << stagger.min_code_length() << "), stagger)\n";
  }
  const auto ctl = bc_D.max_code_length()
                 + 1
                 + bc_B.max_code_length()
                 + bc_A.max_code_length()
                 + 1
                 + bc_C.max_code_length();
  cout << "    code_total_length:  " << ctl << '\n';
  cout << "output\n";
  cout << "    R1:  " << R1_out->path() << '\n';
  cout << "    R2:  " << R2_out->path() << '\n';
}


int main(int argc, const char** argv) {
  try {
    bool force = false;
    bool dry_run = false;
    std::vector<std::string> replace{};
    std::string json_file;
    int i = 1;
    while (i < argc) {
      if (0 == std::strcmp(argv[i], "-h") * std::strcmp(argv[i], "--help")) {
        throw usage_msg;
      }
      else if (0 == std::strcmp(argv[i], "-f") * std::strcmp(argv[i], "--force")) {
        force = true;
      }
      else if (0 == std::strcmp(argv[i], "--dry")) {
        dry_run = true;
      }
      else if (0 == std::strcmp(argv[i], "--replace")) {
        if ((i + 1) > argc) throw "--replace: missing arguments";
        replace.emplace_back(argv[i+1]);
        ++i;
      }
      else if (std::filesystem::exists(argv[i])) {
        json_file = argv[i];
      }
      else {
        std::cerr << "invalid argument '" << argv[i] << "'\n";
        throw usage_msg;
      }
      ++i;
    }
    if (!fs::exists(json_file)) {
      throw std::runtime_error("json file doesn't exists");
    }
    std::ifstream is(json_file);
    auto J = nlohmann::json::parse(is, nullptr, true, /* ignore_comments */ true);
    auto h4 = H4{J};
    if (fs::exists(h4.out_dir)) { 
      if (!force) {
        throw "Output directory already exists. Consider '-f'";
      }
      fs::remove_all(h4.out_dir);
    }
    fs::create_directories(h4.out_dir);
    (dry_run) ? h4.dry_run() : h4.run();

    return 0;
  }
  catch (std::exception& err) {
    std::cerr << err.what() << std::endl;
  }
  catch (const char* err) {
    std::cerr << err << std::endl;
  }
  catch (...) {
    std::cerr << "unknown exception" << std::endl;
  }
  return 1;
}
