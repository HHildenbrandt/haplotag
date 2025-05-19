#pragma once

#include <iostream>
#include <chrono>
#include <functional>
#include <string>


// returns average time per repetition [s]
template <typename Fun>
auto bench(const std::string msg, Fun&& fun, int rep = 1000) {
  auto t0 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < rep; ++i) {
    std::invoke(fun);
  }
  auto time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0);
  std::cout << msg << " passed. " << time.count() << " ms,  " <<  time.count() / rep << " ms per rep";
  return time.count() / (1000.0 * rep);
}
