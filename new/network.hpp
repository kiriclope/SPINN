#ifndef NETWORK_HPP
#define NETWORK_HPP

#include <vector>
#include "globals.hpp"

class RateNetwork {
public:
  std::vector<float> rates;
  std::vector<float> ff_inputs;
  std::vector<float> inputs;
  std::vector<float> net_inputs;
  std::vector<float> dot_prod;
  std::vector<std::vector<float>> matrix;

  RateNetwork() :
    rates(N, 0.0),
    ff_inputs(N, 1.0),
    inputs(N, 0.0),
    net_inputs(N, 0.0),
    dot_prod(N, 0.0),
    matrix(N, std::vector<float>(N, 0.0)) {}

  std::vector<std::vector<float>> generateMatrix();

  void initNetwork();
  void updateInputs();
  void updateRates();
  void runSimul();
};

#endif
