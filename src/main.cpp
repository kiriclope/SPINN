#include <iostream>
#include <iomanip>
#include <chrono>
#include <string>

#include "globals.hpp"
#include "utils.hpp"
#include "sparse_mat.hpp"
#include "lif_network.hpp"

#define YAML_HEADER_ONLY 1

int main(int argc , char** argv) {
  auto start_time = std::chrono::high_resolution_clock::now();

  std::cout << "LIF NETWORK SIMULATION" << std::endl;

  std::string configname = argv[1] ;

  std::cout << "Loading config from: " << configname;
  loadConfig(configname);
  std::cout << " Done" << std::endl;

  init_lif();
  runSimul();
  free_lif();

  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();

  int hours = duration / 3600;
  int minutes = (duration % 3600) / 60;
  int seconds = duration % 60;

  std::cout << "Total time: ";
  std::cout << std::setw(2) << std::setfill('0') << hours << ":"
            << std::setw(2) << std::setfill('0') << minutes << ":"
            << std::setw(2) << std::setfill('0') << seconds << std::endl;

  return 0;
}

