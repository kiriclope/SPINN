#include <cmath>
#include <vector>

#include "network.hpp"
#include "utils.hpp"


std::vector<std::vector<float>> RateNetwork::generateMatrix() {

  std::vector<float> cosines(N/2+1);
  for (int k = 0; k <= N/2; k++) {
    float angle = 2.0 * M_PI * k / N;
    cosines[k] = std::cos(angle);
  }

  std::vector<std::vector<float>> C(N, std::vector<float>(N, 0.0));
  for (int i = 0; i < N; i++) {
    for (int j = i; j < N; j++) {
      int k = std::abs(j-i);
      if (k > N/2) {
        k = N - k;
      }
      C[i][j] = J0 / N * (1.0 + 2.0 * J1 * cosines[k]);
    }
  }
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < i; j++) {
      C[i][j] = C[j][i];
    }
  }
  return C;
}

void RateNetwork::initNetwork(){
  rates = generateGaussianVector<float>(N, 2.0, 0.5);
  inputs = generateGaussianVector<float>(N, 2.0, 0.5);
  transformVector(ff_inputs, I0, std::multiplies<decltype(ff_inputs)::value_type>());
}

void RateNetwork::updateInputs(){
  // transformVector(inputs, EXP_DT_TAU_SYN, std::multiplies<decltype(inputs)::value_type>());

  // dotProduct(matrix, rates, dot_prod);
  // transformVector(dot_prod, DT_TAU_SYN, std::multiplies<decltype(dot_prod)::value_type>());

  // sumVectors(inputs, dot_prod, inputs);

  dotProduct(matrix, rates, inputs);
  sumVectors(ff_inputs, inputs, net_inputs);
}

void RateNetwork::updateRates(){
  transformVector(rates, EXP_DT_TAU_RATE, std::multiplies<decltype(rates)::value_type>());
  transferFunc(net_inputs, THRESH, net_inputs);
  transformVector(net_inputs, DT_TAU_RATE, std::multiplies<decltype(net_inputs)::value_type>());
  sumVectors(net_inputs, rates, rates);
}

void RateNetwork::runSimul(){

  bool append = false;

  matrix = generateMatrix();
  initNetwork();

  for(float time = 0; time < DURATION; time += DT) {
    updateInputs();
    updateRates();

    // if(time % T_WINDOW == 0) {
    //   writeVectorToFile(rates, "rates.txt", append);
    //   writeVectorToFile(inputs, "inputs.txt", append);
    //   if(append == false) append = true;
    // }
  }

}
