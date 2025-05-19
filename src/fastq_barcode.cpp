#include <iostream>
#include <fastq/barcode.hpp>
#include <fastq/fuzzy_matching.hpp>

using namespace fastq;


int main() {
  try {
    auto bc_A = barcode_t("../data/BC_A_H4.txt");
    auto bc_B = barcode_t("../data/BC_B.txt");
    auto bc_C = barcode_t("../data/BC_C.txt");
    auto bc_D = barcode_t("../data/BC_D.txt");


    for (auto& bc : bc_A) {
      std::cout << bc.code_letter() << "  " << bc.tag() << " " << bc.code() << '\n';
    }
    std::cout << "min code length: " << bc_A.min_code_length() << "  max code length: " << bc_A.max_code_length() << '\n';
    return 0;
  }
  catch (std::exception& err) {
    std::cerr << err.what() << std::endl;
  }
  return -1;
}
