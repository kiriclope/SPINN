#include <iostream>
#include <chrono>

#include "globals.hpp"
#include "network.hpp"
#include "utils.hpp"

int main(int argc , char** argv) {
  auto start_time = std::chrono::high_resolution_clock::now();

  RateNetwork net;

  net.initNetwork();

  std::cout << net.rates[0] << std::endl;

  net.runSimul();

  std::cout << net.rates[0] << std::endl;

  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();

  std::cout << "Total time: " << duration << " seconds" << std::endl;

  return 0;
}

