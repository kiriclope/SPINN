#ifndef LIF_STP_HPP
#define LIF_STP_HPP

#include "globals.hpp"

class LifNetwork;

class StpModel {
public:
  float *u_stp;
  float *x_stp;
  float *A_stp;

  float *TAU_FAC;
  float *TAU_REC;
  float USE;

  StpModel():
    u_stp(new float[N]()),
    x_stp(new float[N]()),
    A_stp(new float[N]()),

    TAU_FAC(new float[N_POP]()),
    TAU_REC(new float[N_POP]())
  {}

  // Destructor to free memory
  ~StpModel() {
    delete[] u_stp;
    delete[] x_stp;
    delete[] A_stp;
  }

  void matoModel(int i, float ISI);
};

#endif
