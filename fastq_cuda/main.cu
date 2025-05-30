#include <iostream>
#include <type_traits>
#include <cuda_runtime.h>
#include <cuda/std/span>
#include <fastq/fastq.hpp>
#include <fastq/reader.hpp>
#include <fastq/splitter.hpp>
#include <fastq/barcode.hpp>


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
      std::cout << "0x" << std::hex << uintptr_t(ptr) << std::dec << " + " << bytes << std::endl;
      cudaDeviceSynchronize();
      // cudaMemAdvise(ptr, bytes, cudaMemAdviseSetPreferredLocation, cudaCpuDeviceId);
      // cudaMemAdvise(ptr, bytes, cudaMemAdviseSetAccessedBy, cudaCpuDeviceId);
      return ptr;
    };

    static void free(void* ptr) noexcept { 
      std::cout << "0x" << std::hex << uintptr_t(ptr) << std::dec << " - " << std::endl;
      cudaDeviceSynchronize();
      cudaFree(ptr); 
    } 
  };


  using reader_t = fastq::detail::reader_t<chunk_allocator<char>, 16 * 1024, 16 * 1024 * 1024>;


  __host__ __device__
  struct str_view {
    const char* first;
    size_t length;
  };


  struct read_line_split {
    using value_type = str_view;

    static value_type apply(fastq::str_view& /* in out */ cv) noexcept {
      using field_split = fastq::policy::delim_split<char, '\n', 0, 1>;
      auto ret = value_type{};
      for (auto i = 0; i < 4; ++i) {
        auto field = field_split::apply(cv);
        if (i == 1) {
          ret = { field.begin(), field.length() };
        }
      }
      return ret;
    }
  };


  using read_field_splitter = fastq::base_splitter<
    cu::reader_t,
    fastq::chunk_splitter<
      fastq::policy::delim_chunk_trim<const char*, fastq::SeqDelim>,
      read_line_split
    >
  >;


  struct barcode_split {
    using value_type = str_view;

    static value_type apply(auto /* in out */ cv) noexcept {
      using line_split = fastq::policy::delim_split<char, '\n', 0, 1>;
      auto line = line_split::apply(cv);
      auto code = line.substr(line.find_last_of(" \t") + 1);
      return { code.begin(), code.length() };
    }
  };


  using barcode_splitter = fastq::base_splitter<
    reader_t,
    fastq::chunk_splitter<
      fastq::policy::delim_chunk_trim<char, '\n'>,
      barcode_split
    >
  >;


  template <typename T>
  auto make_device_ptr(const T* hptr, size_t n) {
    void* dptr = nullptr;
    if (auto res = cudaMalloc(&dptr, n * sizeof(T)); cudaSuccess != res) {
      throw std::bad_alloc();
    }
    cudaMemcpy(dptr, hptr, n * sizeof(T), cudaMemcpyHostToDevice);
    return std::shared_ptr<T[]>((T*)dptr, cudaFree);
  }


  template <size_t N>
  __device__
  size_t edit_distance(const char* a, size_t m, const char* b, size_t n) {
    if (m > n) {
      auto s = b; b = a; a = s; 
      auto l = n; n = m; m = l; 
    }
    // remove matching prefixes and suffixes
    while (m && (*a == *b)) { ++a; ++b; --m; --n; }
    while (m && (a[m-1] == b[n-1])) { --m; --n; }
    size_t D[N + 1];  // scratch
    for (size_t i = 0; i < m + 1; ++i) D[i] = i;
    for (auto i = 1; i <= n; ++i) {
      const auto bi = b[i - 1];
      auto tmp = D[0]; D[0] = i;
      for (auto j = 1; j <= m; ++j) {
        if (a[j - 1] != bi) {
          tmp = min(D[j], min(D[j - 1], tmp)) + 1;
        }
        auto t = tmp; tmp = D[0]; D[0] = t;
      }
    }
    return D[m];
  }


  __global__
  void blk_edit_distance(str_view* RX, size_t nrx, str_view* BC, size_t nbc, size_t* out) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx > nrx) return;
    str_view rx = RX[idx];
    for (int i = 0; i < nbc; ++i) {
      str_view bc = BC[i];
      out[idx] = edit_distance<8>(rx.first, 7, bc.first, bc.length);
    }
  }

  __global__
  void nop(str_view* RX, size_t nrx, str_view* BC, size_t nbc, size_t* out) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx > nrx) return;
    out[idx] = idx;
  }

  void launch_blk_edit_distance(auto& blk, str_view* dbc, size_t nbc, size_t* ed_out) {
    auto drx = make_device_ptr<str_view>(blk.data(), blk.size());
    dim3 numBlocks(blk.size() / 256);
    blk_edit_distance<<<numBlocks, 256>>>(drx.get(), blk.size(), dbc, nbc, ed_out);
  }

}


int main() {
  constexpr size_t N = 10000;
  auto s = cu::read_field_splitter("../data/_gen_I2_001.fastq.gz");
  auto bc = cu::barcode_splitter("../data/BC_A.txt")(N);
  auto dbc = cu::make_device_ptr(bc.data(), bc.size());
  void* ed_out = nullptr;
  cudaMalloc(&ed_out, N * sizeof(cu::str_view));
  size_t items = 0;
  while (!s.eof()) {
    auto blk = s(N);
    cu::launch_blk_edit_distance(blk, dbc.get(), bc.size(), (size_t*)ed_out);
    items += blk.size();
  }
  std::cout << items << " items read, " << s.reader().tot_bytes() / (1000 * 1000) << " MB" << std::endl;
  return 0;
}
