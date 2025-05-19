#include <type_traits>
#include <cuda_runtime.h>
#include <cuda/std/span>
#include <fastq/reader.hpp>
#include <fastq/splitter.hpp>


namespace cu {

  template <typename T>
  struct chunk_allocator {
    static_assert(std::is_trivially_constructible_v<T>);
    static_assert(std::is_trivially_destructible_v<T>);

    static void* alloc(size_t bytes) { 
      void* ptr; 
      if (cudaSuccess != cudaMallocManaged(&ptr, bytes)) {
        throw std::bad_alloc{};
      }
      cudaMemAdvise(ptr, bytes, cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId);
      cudaMemAdvise(ptr, bytes, cudaMemAdviseSetAccessedBy, cudaCpuDeviceId);
      return ptr;
    };

    static void free(void* ptr) noexcept { cudaFree(ptr); } 
  };


  using reader_t = fastq::detail::reader_t<chunk_allocator<char>, 64 * 1024 * 1024>;
  //using reader_t = fastq::reader_t;

}

using splitter_t = fastq::seq_splitter<cu::reader_t>;


int main() {
  auto s = splitter_t("../data/_gen_I2_001.fastq.gz");
  while (!s.eof()) {
    auto x = s();
  }
  return 0;
}
