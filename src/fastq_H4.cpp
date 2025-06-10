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
      throw std::runtime_error("json files doesn't exists");
    }
    auto json_dir = json_file;
    json_dir.remove_filename();
    std::ifstream is(json_file);
    auto J = json::parse(is, nullptr, true, /* ignore_comments */ true);
    auto jin = J.at("input");
    data_dir = json_dir / jin.at("dir").get<fs::path>();
    auto jbc = jin.at("barcodes");
    auto gen_bc = [&](const char* L) {
      return fastq::barcode_t{data_dir / jbc.at(L).at("file").get<fs::path>(), jbc.at(L).at("unclear_tag")};
    };
    bc_A = gen_bc("A");
    bc_B = gen_bc("B");
    bc_C = gen_bc("C");
    bc_D = gen_bc("D");
    optional_json(plate = gen_bc("plate"));
    optional_json(stagger = gen_bc("stagger"));
    R1 = Rsplitter{data_dir / jin.at("R1").get<fs::path>()};
    optional_json(R2 = Rsplitter{data_dir / jin.at("R2").get<fs::path>()});
    R3 = Rsplitter{data_dir / jin.at("R3").get<fs::path>()};
    optional_json(I1 = Isplitter{data_dir / jin.at("I1").get<fs::path>()});
    I2 = Isplitter{data_dir / jin.at("I2").get<fs::path>()};
    
    // sanity checks
    if (1 == (I1.failed() + plate.empty())) {
      throw std::runtime_error("only one of 'plate', 'I1' given");
    }
    if (1 == (R2.failed() + stagger.empty())) {
      throw std::runtime_error("only one of 'stagger', 'R2' given");
    }
    
    auto jout = J.at("output");    
    out_dir = json_dir / jout.at("dir").get<fs::path>();
    R1_out.reset(new fastq::writer_t<>{out_dir/ jout.at("R1").get<fs::path>(), gPool});
    R2_out.reset(new fastq::writer_t<>{out_dir/ jout.at("R2").get<fs::path>(), gPool});
  }

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
  fs::path data_dir;
  fs::path out_dir;
};


void dry_run(H4 h4) {
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
  cout << "fastq.gz files:\n";
  gz_stats("    R1", h4.R1);
  gz_stats("    R2", h4.R2);
  gz_stats("    R3", h4.R3);
  gz_stats("    I1", h4.I1);
  gz_stats("    I2", h4.I2);

  cout << "matches\n";
  const auto ctl = h4.bc_D.max_code_length()
                 + 1
                 + h4.bc_B.max_code_length()
                 + h4.bc_A.max_code_length()
                 + 1
                 + h4.bc_C.max_code_length();
  cout << "  code_total_length: " << ctl << '\n';
}


int main(int argc, const char** argv) {
  try {
    fs::path json_file{"../src/H4.json"}; // default for debugging
    if (argc > 1) {
      json_file = argv[1];
    }
    dry_run(H4(fs::current_path() / json_file));
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
