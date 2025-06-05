#include <iostream>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <fastq/barcode.hpp>
#include <fastq/reader.hpp>
#include <fastq/splitter.hpp>


struct barcode_json {
  std::filesystem::path file;
  std::string unclear_tag;
  auto value(const std::filesystem::path& parent) { return fastq::barcode_t(parent / file, unclear_tag); }
};

void from_json(const nlohmann::json& J, barcode_json& val) {
  val.file = J.at("file").get<std::filesystem::path>();
  val.unclear_tag = J.at("unclear_tag").get<std::string>();
}


struct H4 {
  using Rsplitter = fastq::line_splitter<>;
  using Isplitter = fastq::seq_field_splitter<0b1010>;   // [RX, QX]

  explicit H4(nlohmann::json& J) {
    data_dir = J.at("data_dir").get<std::filesystem::path>();
    auto bcJ = J.at("barcodes");
    bc_A = bcJ.at("A").get<barcode_json>().value(data_dir);
    bc_B = bcJ.at("B").get<barcode_json>().value(data_dir);
    bc_C = bcJ.at("C").get<barcode_json>().value(data_dir);
    bc_D = bcJ.at("D").get<barcode_json>().value(data_dir);
    try { // optional
      plate = bcJ.at("plate").get<barcode_json>().value(data_dir);
    }
    catch (...) {}
    stagger = bcJ.at("stagger").get<barcode_json>().value(data_dir);
    auto ji = J.at("in");
    R1 = std::make_unique<Rsplitter>(data_dir / ji.at("R1").get<std::filesystem::path>());
    R2 = std::make_unique<Rsplitter>(data_dir / ji.at("R2").get<std::filesystem::path>());
    R3 = std::make_unique<Rsplitter>(data_dir / ji.at("R3").get<std::filesystem::path>());
    I1 = std::make_unique<Isplitter>(data_dir / ji.at("I1").get<std::filesystem::path>());
    I2 = std::make_unique<Isplitter>(data_dir / ji.at("I2").get<std::filesystem::path>());
  }

  std::filesystem::path data_dir;
  fastq::barcode_t bc_A;
  fastq::barcode_t bc_B;
  fastq::barcode_t bc_C;
  fastq::barcode_t bc_D;
  fastq::barcode_t plate;
  fastq::barcode_t stagger;

  std::unique_ptr<Rsplitter> R1;
  std::unique_ptr<Rsplitter> R2;
  std::unique_ptr<Rsplitter> R3;
  std::unique_ptr<Isplitter> I1;
  std::unique_ptr<Isplitter> I2;
};

 

int main(int argc, const char** argv) {
  try {
    std::filesystem::path json_file{"../src/H4.json"}; // default for debugging
    if (argc > 1) {
      json_file = argv[1];
    }
    if (!std::filesystem::exists(json_file)) {
      throw std::runtime_error("json file doesn't exists");
    }
    std::ifstream is(json_file);
    auto J = nlohmann::json{};
    is >> J;
    auto h4 = H4(J);
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