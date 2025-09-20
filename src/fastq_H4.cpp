#include <iostream>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <fastq/barcode.hpp>
#include <fastq/reader.hpp>
#include <fastq/splitter.hpp>
#include <fastq/writer.hpp>
#include "device/pool.hpp"


namespace fs = std::filesystem;
using json = nlohmann::json;


// global thread-pool
std::shared_ptr<hahi::pool_t> gPool{ new hahi::pool_t{unsigned(-1)} };


#define optional_json(expr) try { expr; } catch (json::exception&) {}


struct H4 {
  using Rsplitter = fastq::line_splitter<>;
  using Isplitter = fastq::seq_field_splitter<0b1010>;   // [RX, QX]

  explicit H4(const fs::path& json_file) {
    if (!fs::exists(json_file)) {
      throw std::runtime_error("json file doesn't exists");
    }
    std::ifstream is(json_file);
    auto J = json::parse(is, nullptr, true, /* ignore_comments */ true);
    auto root = J.at("root").get<fs::path>();
    std::cout << root << '\n';
    auto jin = J.at("input");
    auto jbc = jin.at("barcodes");
    auto bc_prefix = root / jbc.at("prefix").get<std::string>();
    std::cout << bc_prefix << '\n';
    auto gen_bc = [&](const char* L) {
      auto bc = fastq::barcode_t{bc_prefix / jbc.at(L).at("file").get<std::string>(), jbc.at(L).at("unclear_tag")};
      optional_json(bc.reset_code_letter(jbc.at(L).at("code_letter").get<std::string>()[0]));
      return bc;
    };
    bc_A = gen_bc("A");
    bc_B = gen_bc("B");
    bc_C = gen_bc("C");
    bc_D = gen_bc("D");
    optional_json(plate = gen_bc("plate"));
    optional_json(stagger = gen_bc("stagger"));
    auto jgz = jin.at("fastq");
    gz_prefix = root.string() + jgz.at("prefix").get<std::string>();
    R1 = Rsplitter{gz_prefix + jgz.at("R1").get<std::string>()};
    optional_json(R2 = Rsplitter{gz_prefix + jgz.at("R2").get<std::string>()});
    R3 = Rsplitter{gz_prefix + jgz.at("R3").get<std::string>()};
    optional_json(I1 = Isplitter{gz_prefix + jgz.at("I1").get<std::string>()});
    I2 = Isplitter{gz_prefix + jgz.at("I2").get<std::string>()};
    
    // sanity checks
    if (1 == (I1.failed() + plate.empty())) {
      throw std::runtime_error("only one of 'plate', 'I1' given");
    }
    if (1 == (R2.failed() + stagger.empty())) {
      throw std::runtime_error("only one of 'stagger', 'R2' given");
    }
   
    auto jout = J.at("output");    
    out_prefix = jout.at("prefix").get<std::string>();
    R1_out.reset(new fastq::writer_t<>{out_prefix + jout.at("R1").get<std::string>(), gPool});
    R2_out.reset(new fastq::writer_t<>{out_prefix + jout.at("R2").get<std::string>(), gPool});
  }
 
  bool has_stagger() const noexcept { return !stagger.empty(); }
  bool has_plate() const noexcept { return !plate.empty(); }

  fastq::barcode_t bc_A;
  fastq::barcode_t bc_B;
  fastq::barcode_t bc_C;
  fastq::barcode_t bc_D;
  fastq::barcode_t plate;
  fastq::barcode_t stagger;

  Rsplitter R1;
  Rsplitter R2;
  Rsplitter R3;
  Isplitter I1;
  Isplitter I2;

  std::unique_ptr<fastq::writer_t<>> R1_out;
  std::unique_ptr<fastq::writer_t<>> R2_out;

  bool clipping = false;
  std::string bc_prefix;
  std::string gz_prefix;
  std::string out_prefix;
};


void dry_run(const H4& h4) {
  using std::cout;
  auto bc_stats = [](const char* name, const auto& bc) { 
    cout << name << "  ";
    if (bc.empty()) {
      cout << "NA\n";
      return;
    }
    cout << '"' << bc[0].tag << "\"  "
         << bc.size() << "  "
         << '[' << bc.min_code_length() << ", " << bc.max_code_length() << "]  "
         << bc.path()
         << '\n';
  };
  cout << "barcodes\n";
  bc_stats("    bc_A   ", h4.bc_A);
  bc_stats("    bc_B   ", h4.bc_B);
  bc_stats("    bc_C   ", h4.bc_C);
  bc_stats("    bc_D   ", h4.bc_D);
  bc_stats("    plate  ", h4.plate);
  bc_stats("    stagger", h4.stagger);

  auto gz_stats = [](const char* name, const auto& gz) {
    cout << name << "  ";
    if (gz.failed()) cout << "NA\n";
    else cout << gz.reader().path() << '\n';
  };
  cout << "fastq.gz files\n";
  gz_stats("    R1", h4.R1);
  gz_stats("    R2", h4.R2);
  gz_stats("    R3", h4.R3);
  gz_stats("    I1", h4.I1);
  gz_stats("    I2", h4.I2);

  cout << "matches\n";
  if (h4.has_stagger()) {
    std::cout << "    S <- idx min_ed(R2[1](0:" << h4.stagger.min_code_length() << "), stagger)\n";
  }
  const auto ctl = h4.bc_D.max_code_length()
                 + 1
                 + h4.bc_B.max_code_length()
                 + h4.bc_A.max_code_length()
                 + 1
                 + h4.bc_C.max_code_length();
  cout << "    code_total_length  " << ctl << '\n';
  cout << "output\n";
  cout << "    R1  " << h4.R1_out->path() << '\n';
  cout << "    R2  " << h4.R2_out->path() << '\n';
}


int main(int argc, const char** argv) {
  try {
    fs::path json_file{"../src/H4.json"}; // default for debugging
    if (argc > 1) {
      json_file = argv[1];
    }
    auto h4 = H4{fs::current_path() / json_file};
    dry_run(h4);
    return 0;
  }
  catch (std::exception& err) {
    std::cerr << err.what() << std::endl;
  }
  catch (...) {
    std::cerr << "unknown exception" << std::endl;
  }
  return 1;
}
