#include <cstdlib>
#include <charconv>
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


constexpr char usage_msg[] = R"(Usage: fastq_h4 JSON_FILE [OPTIONS]...
  -h, --help: show this message.
  -f, --force: force overwrite of output directory.
  -v, --verbose: verbose output.
  --replace '{"json_pointer": value}'.
    Ex: --replace '{"/range": "0-1000"}' --replace '{"/barcode/plate/file": "Plate_BC_7.txt"}'
  --dry: dry-run.
)";


namespace fs = std::filesystem;
using namespace std::chrono_literals;
using json = nlohmann::json;
using Splitter = fastq::seq_field_splitter<0b1111>;


// global thread-pool
// the pool is shared betweeen matching and compressing jobs
std::shared_ptr<hahi::pool_t> gPool;


// hacky way replaces leadiing "~/" in path with the path to HOME 
fs::path expand_home(const fs::path& path) {
  // this is a bit hacky
  if (path.string().starts_with("~/")) {
    char* home = nullptr;
#ifdef _WIN32_      
    home = std::getenv("HOMEPATH");
#else
    home = std::getenv("HOME");
#endif
    if (home) return std::filesystem::path(home) / path.string().substr(2);
  }
  return path;
}


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


// give reads a name...
enum ReadIdx {
  R1_, R2_, R3_, R4_, I1_
};


#define optional_json(expr) try { expr; } catch (json::exception&) {}


struct H4 {
  using blks_t = std::vector<Splitter::blk_type>; 
  static constexpr size_t blk_size = 10000;   // non-sensitive tuneable

  explicit H4(const json& Jin, bool verbose) : verbose(verbose), J(Jin) {
    range = parse_range(J.at("range").get<std::string>());
    if (range.first >= range.second) throw "invalid range";

    // create thread pool
    gPool.reset( new hahi::pool_t(120/*J.at("pool_threads").get<unsigned>()*/));

    // barcodes
    auto jbc = J.at("barcodes");
    auto bc_root = expand_home(jbc.at("root").get<fs::path>());
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
    // plate is optional
    if (!jbc.at("plate").at("file").get<std::string>().empty()) {
      plate = gen_bc("plate");
    }
    stagger = gen_bc("stagger");

    // reads
    auto jr = J.at("reads");
    gz_root = expand_home(jr.at("root").get<std::string>());
    R1 = Splitter{gz_root / jr.at("R1").get<std::string>()};
    R2 = Splitter{gz_root / jr.at("R2").get<std::string>()};
    R3 = Splitter{gz_root / jr.at("R3").get<std::string>()};
    R4 = Splitter{gz_root / jr.at("R4").get<std::string>()};
    if (!plate.empty()) {
      I1 = Splitter{gz_root / jr.at("I1").get<std::string>()};
    }
    // output
    auto jout = J.at("output"); 
    r1_out = !jout.at("R1").get<std::string>().empty();  
    clipping = !jout.at("R2").get<std::string>().empty();
    out_root = expand_home(jout.at("root").get<std::string>());
  }

  bool has_stagger() const noexcept { return !stagger.empty(); }
  bool has_plate() const noexcept { return !plate.empty(); }

  void dry_run() {
    using std::cout;
    cout << "range: " << range.first << '-' << range.second << '\n';
    cout << "pool_threads: " << gPool->num_threads() << '\n';
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
    cout << "    R1: " << (r1_out ? out_root / J.at("/output/R1"_json_pointer).get<std::string>() : "NA") << '\n';
    cout << "    R2: " << (clipping ? out_root / J.at("/output/R2"_json_pointer).get<std::string>() : "NA (no clipping)") << '\n';
  }

  // feeding the work pipeline
  template <bool has_plate>
  void run() {
    // layzy creation of writers
    if (r1_out) R1_out.reset(new fastq::writer_t<>{out_root / J.at("/output/R1"_json_pointer).get<std::string>(), gPool});
    if (clipping) R2_out.reset(new fastq::writer_t<>{out_root / J.at("/output/R2"_json_pointer).get<std::string>(), gPool});
    
    std::vector<Splitter*> RS = { &R1, &R2, &R3, &R4 };
    if constexpr (has_plate) RS.push_back(&I1);

    size_t i = 0; // sequence number
    // skip head of range
    for (; i < range.first; ++i) {
      for (auto* R : RS) {
        if (R->eof()) break;
        R->operator()();
      }
    }
    if (i != range.first) throw "range exceeds number of reads";

    auto match_queue = std::deque<std::future<h4_matches_t>>{};
    bool any_eof = false;
    for (; !any_eof && (i < range.second); i += blk_size) {
      // collect block of reads
      auto blks = blks_t{};
      const auto n = std::min(range.second - i, blk_size);  // sequences to read
      for (auto* R : RS) {
        blks.emplace_back(R->operator()(n));
        any_eof |= R->eof();
      }
      // check if all read files had equal sequences
      const auto exp = blks[0].size();
      for (size_t i = 1; i < blks.size() - !has_plate; ++i) {
        if (blks[i].size() != exp) throw "inconsistent number of sequences in input";
      }
      // enqueue block-matching job to one of the matching-thread
      // blocks until a thread is available in gPool
      match_queue.emplace_back(gPool->async([this, blks = std::move(blks)]() mutable {
        return this->blk_match<has_plate>(std::move(blks));
      }));
      // anything ready yet?
      // never blocks
      while (!match_queue.empty() && (std::future_status::ready == match_queue.front().wait_for(0s))) {
        auto m = match_queue.front().get();   // doesn's block either
        match_queue.pop_front();
        process_matches<has_plate>(m);
      }
    }
    // left-overs
    while (!match_queue.empty()) {
      auto m = match_queue.front().get();   // might block
      match_queue.pop_front();
      process_matches<has_plate>(m);
    }
    // dump json to output folder for reference
    auto js = std::ofstream(out_root / "H4.json");
    js << J.dump();
  }

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
  bool r1_out = false;
  std::filesystem::path bc_root;
  std::filesystem::path gz_root;
  std::filesystem::path out_root;
  const json& J;

private:
  struct h4_match_t {
    int sn = 0;
    fastq::match_t s, a, b, c, d, p;
    bool any_invalid = false;
    bool any_unclear = false;
  };
  using h4_matches_t = std::pair<std::vector<h4_match_t>, blks_t>;

  // matching 
  template <bool has_plate>
  h4_matches_t blk_match(blks_t&& blks) {
    auto matches = std::vector<h4_match_t>{};
    matches.reserve(blks.size());
    // expected code lengths
    const size_t scl = stagger.max_code_length();
    const size_t bcl = bc_B.max_code_length();           
    const size_t dcl = bc_D.max_code_length();
    const size_t ccl = bc_C.max_code_length();  
    const size_t pcl = plate.max_code_length();   // always defined, 0 if empty
    std::string RX{};
    for (size_t i = 0; i < blks[0].size(); ++i) {
      auto& m = matches.emplace_back();
      m.s = fastq::min_edit_distance(fastq::max_substr(blks[R4_][i][1], 0, scl), scl, stagger);
      m.sn = (m.s.rt <= fastq::ReadType::unclear) ? 0 : m.s.idx - 1;   // rerquires 'sorted' stagger barcodes
      RX = blks[R2_][i][1];
      RX.append(blks[R3_][i][1]);
      m.b = fastq::min_edit_distance(fastq::max_substr(RX, bcl + 1, bcl), bcl, bc_B);
      m.d = fastq::min_edit_distance(fastq::max_substr(RX, 0, dcl), dcl, bc_D);
      const auto acl = bc_A.min_code_length() + m.sn;
      m.a = fastq::min_edit_distance(fastq::max_substr(RX, bcl + dcl + 1, acl), acl, bc_A);
      m.c = fastq::min_edit_distance(fastq::max_substr(RX, bcl + dcl + acl + 2, ccl), ccl, bc_C);
      if constexpr (has_plate) {
        m.p = fastq::min_edit_distance(fastq::max_substr(blks[I1_][i][1], 0, pcl), pcl, plate);
        m.any_invalid = (m.p.rt == fastq::ReadType::invalid);
        m.any_unclear = (m.p.rt == fastq::ReadType::unclear);
      }
      // summary
      for (fastq::ReadType rt : { m.s.rt, m.a.rt, m.b.rt, m.c.rt, m.d.rt }) {
        m.any_invalid |= (rt == fastq::ReadType::invalid);
        m.any_unclear |= (rt == fastq::ReadType::unclear);
      }
    }
    return { std::move(matches), std::move(blks) };
  }

  // output processing
  template <bool has_plate, bool has_clipping>
  void do_process_matches(const h4_matches_t& h4_matches,
                          auto put,
                          auto puts) {
    const auto& matches = h4_matches.first;
    const auto& blks = h4_matches.second;
    const size_t pcl = plate.max_code_length();   // always defined, 0 if empty
    for (size_t i = 0; i < blks[0].size(); ++i) {
      const auto& match = matches[i];

      // compile comments
      const auto name = blks[R1_][i][0];
      put(name.substr(0, name.find_first_of(" \t")));
      put("\tBX:Z:");
      put(bc_A[match.a.idx].tag);
      put(bc_C[match.c.idx].tag);
      put(bc_B[match.b.idx].tag);
      put(bc_D[match.d.idx].tag);
      if constexpr (has_plate) {
        put("-");
        put(plate[match.p.idx].tag);
      }
      put("\tRX:Z:");
      put(blks[R2_][i][1]);
      put(blks[R3_][i][1]);
      if constexpr (has_plate) {
        put("+");
        put(blks[I1_][i][1]);
      }
      put("\tQX:Z:");
      put(blks[R2_][i][3]);
      put(blks[R3_][i][3]);
      if constexpr (has_plate) {
        put("+");
        put(blks[I1_][i][3]);
      }
      puts({});

      // copy over unchanged fields to R1_out
      for (auto j : {1,2,3}) R1_out->puts(blks[R1_][i][j]);

      if constexpr (has_clipping) {
        // copy clipped fields to R2_out
        auto clip_size = stagger.max_code_length() + 1;
        clip_size += (match.a.rt == fastq::ReadType::unclear) 
                    ? bc_A.max_code_length()
                    : bc_A[match.a.idx].code.length();
        R2_out->puts(fastq::max_substr(blks[R4_][i][1], clip_size));
        R2_out->puts(blks[R4_][i][2]);
        R2_out->puts(fastq::max_substr(blks[R4_][i][3], clip_size));
      }
    }
  }

  template <bool has_plate>
  void process_matches(const h4_matches_t& h4_matches) {
    if (clipping) {
      // creates a lot of redundant data...
      do_process_matches<has_plate, true>(
        h4_matches, 
        [this](fastq::str_view str) { R1_out->put(str); R2_out->put(str); },
        [this](fastq::str_view str) { R1_out->puts(str); R2_out->puts(str); }
      );
    }                                   
    else {
      do_process_matches<has_plate, false>(
        h4_matches, 
        [this](fastq::str_view str) { R1_out->put(str); },
        [this](fastq::str_view str) { R1_out->puts(str); }
      );
    }
  }
};


int main(int argc, const char** argv) {
  try {
    bool force = false;
    bool verbose = false;
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
      else if (0 == std::strcmp(argv[i], "-v") * std::strcmp(argv[i], "--verbose")) {
        verbose = true;
      }
      else if (0 == std::strcmp(argv[i], "--dry")) {
        dry_run = true;
      }
      else if (0 == std::strcmp(argv[i], "--replace")) {
        if ((i + 1) > argc) throw "--replace: missing arguments";
        replace.emplace_back(argv[++i]);
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
    for (auto& r : replace) {
      auto R = json::parse(r);
      for (auto& e : R.items()) {
        J[json::json_pointer(e.key())] = e.value();
      }
    }
    auto h4 = H4{J, verbose};
    if (dry_run) {
      h4.dry_run();
      return 0;
    }
    if (!(h4.clipping || h4.r1_out)) {
      throw ("Neigther R1 nor R2 output specified\n. Bailing out.");
    }
    if (fs::exists(h4.out_root)) { 
      if (!force) {
        throw "Output directory already exists. Consider '-f'";
      }
      fs::remove_all(h4.out_root);
    }
    fs::create_directories(h4.out_root);
    if (h4.plate.empty()) {
      h4.run<false>();
    }
    else {
      h4.run<true>();
    }
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
