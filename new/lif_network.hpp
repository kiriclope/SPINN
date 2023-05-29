#ifndef NETWORK_HPP
#define NETWORK_HPP
#include "globals.hpp"

class LifNetwork {
public:
  float *rates;
  float *volts;
  float *spikes;
  float *spike_times;

  float *ff_inputs;
  float **inputs;
  float *net_inputs;

  float *Jab_scaled;
  float *Iext_scaled;

  unsigned long *colptr;
  int *indices;

  float *x_stp;
  float *u_stp;
  float *A_stp;

  LifNetwork(): rates(new float[N]()),
                volts(new float[N]()),
                spikes(new float[N]()),
                spike_times(new float[N]()),

                ff_inputs(new float[N]()),
                inputs(new float*[N_POP]()),
                net_inputs(new float[N]()),

                Jab_scaled(new float[N_POP * N_POP]()),
                Iext_scaled(new float[N_POP]()),

                colptr(new unsigned long[N+1]()),
                indices(new int[(unsigned long) (N * 5.0 * K)]()),

                x_stp(new float[N]()),
                u_stp(new float[N]()),
                A_stp(new float[N]())
  {
    for (int i = 0; i < N_POP; i++) {
      inputs[i] = new float[N]();
    }

    if(IF_STP) {
      x_stp = new float[N]();
      u_stp = new float[N]();
      A_stp = new float[N]();
    }
  }

  // Destructor to free memory
  ~LifNetwork() {
    delete[] rates;
    delete[] volts;
    delete[] spikes;
    delete[] spike_times;

    delete[] ff_inputs;

    for (int i = 0; i < N_POP; i++) {
      delete[] inputs[i];
    }
    delete[] inputs;

    delete[] net_inputs;
    delete[] Jab_scaled;
    delete[] Iext_scaled;

    delete[] colptr;
    delete[] indices;

    if (IF_STP) {
      delete[] x_stp;
      delete[] u_stp;
      delete[] A_stp;
    }
  }

  void printParam();
  void initNetwork();
  void updateVolts();
  void updateFFinputs(int step);
  void updateInputs();
  void updateSpikes(int step);
  void updateStp(int i, int step);
  void runSimul();
};

#endif
