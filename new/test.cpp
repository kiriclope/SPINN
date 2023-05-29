#include<iostream>
#include <vector>
#include <yaml-cpp/yaml.h>

#include "test.hpp"
#include "other_test.hpp"
#include "const.hpp"

YAML::Node config;
// int N;
// std::vector<float> Jab;

int main() {
  config = YAML::LoadFile("config_EI.yml");

  N = config["N"].as<int>();
  Jab = config["Jab"].as<std::vector<float>>();

  std::cout << "N " << N << std::endl;
  std::cout << "Jab " << Jab[0]  << " " << Jab[1] << std::endl;

  print_params();

  return 0;
}
