#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>
#include <chrono>

#include "globals.hpp"
#include "utils.hpp"
#include "sparse_mat.hpp"
#include "lif_network.hpp"

float *rates;
float *volts;
float *spikes;
float *spike_times;

float *ff_inputs;
float **inputs;
float *inputs_NMDA;
float *net_inputs;

float *Jab_scaled;
float *Jab_NMDA;

float *Iext_scaled;

size_t *colptr;
int *indices;

float *x_stp;
float *u_stp;
float *A_stp;

std::mt19937 gen;
std::normal_distribution<float> white(0.0, 1.0);

void init_lif() {
  rates = new float[N]();
  volts = new float[N]();
  spikes = new float[N]();
  spike_times= new float[N]();

  ff_inputs = new float[N]();
  inputs = new float*[N_POP]();
  net_inputs = new float[N]();

  Jab_scaled = new float[N_POP * N_POP]();
  Iext_scaled = new float[N_POP]();

  colptr = new size_t[N+1]();
  indices = new int[(size_t) (N * 5.0 * K)]();

  for (int i = 0; i < N_POP; i++)
    inputs[i] = new float[N]();

  if (IF_NMDA) {
    Jab_NMDA = new float[N_POP]();
    inputs_NMDA = new float[N]();
  }

  if(IF_STP) {
    x_stp = new float[N] {1.0};
      u_stp = new float[N] {USE[0]};
      A_stp = new float[N] {USE[0]};
  }

  ensureDirExists(DATA_PATH);
}

// Destructor to free memory
void free_lif() {
  delete[] rates;
  delete[] volts;
  delete[] spikes;
  delete[] spike_times;

  delete[] ff_inputs;

  for (int i = 0; i < N_POP; i++) {
    delete[] inputs[i];
  }
  delete[] inputs;

  if (IF_NMDA) {
    delete[] Jab_NMDA;
    delete[] inputs_NMDA;
  }

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

void initNetwork() {
  // Here (V_THRESH - V_REST) is important to get rates from the balanced linear eqs m = I/J;
  std::cout << "Initializing Network, output dir:" << DATA_PATH;

  // volts = generateGaussianVector<float>(N, 2.0, 0.5);
  for(int i=0; i < N_POP; i++)
    for(int j=0; j < N_POP; j++)
      // Jab_scaled[j + i * N_POP] = GAIN * Jab[j + i * N_POP];
      // Jab_scaled[j + i * N_POP] = GAIN * Jab[j + i * N_POP] / TAU_SYN[j] / Ka[j] * sqrt(Ka[0]);
      Jab_scaled[j + i * N_POP] = GAIN * Jab[j + i * N_POP] * (V_THRESH - V_REST) / TAU_SYN[j] / sqrt(Ka[j]);
      // Jab_scaled[j + i * N_POP] = GAIN * Jab[j + i * N_POP] * (V_THRESH - V_REST) / TAU_SYN[j] / sqrt(K);

  if(IF_NMDA)
    for(int i=0; i < N_POP; i++)
      // Jab_NMDA[i] = GAIN * Jab[i * N_POP] ;
      // Jab_NMDA[i] = GAIN * Jab[i * N_POP] / TAU_NMDA[i] / Ka[0] * sqrt(Ka[0]);
      Jab_NMDA[i] = GAIN * Jab[i * N_POP] * (V_THRESH - V_REST) / TAU_NMDA[i] / sqrt(Ka[0]);

  Jab_NMDA[0] *= (1.0 - R_NMDA[0]) / R_NMDA[0];
  Jab_NMDA[1] *= (1.0 - R_NMDA[1]) / R_NMDA[1];

  for(int i=0; i < N_POP; i++)
    // Iext_scaled[i] = GAIN * Iext[i] ;
    // // Iext_scaled[i] = GAIN * Iext[i] * sqrt(Ka[0]) ;
    Iext_scaled[i] = GAIN * Iext[i] * sqrt(Ka[0]) * (V_THRESH - V_REST);

  for(int i=0; i<N; i++)
    ff_inputs[i] = Iext_scaled[which_pop[i]];

  std::cout << " Done" << std::endl;
}

void updateFFinputs(int step) {

  float theta_i = 0;

  if (step == (int) (T_STIM[0] / DT)) {
    if (VERBOSE)
      std::cout << " STIM ON" << std::endl;
    for (int i = 0; i < N; i++) {
      theta_i = (2.0 * M_PI * (i - cNa[which_pop[i]])) / (float) Na[which_pop[i]];

      ff_inputs[i] = Iext_scaled[which_pop[i]]
        + A_STIM[which_pop[i]] * sqrt(Ka[0])
        * (1.0 + KAPPA_STIM[which_pop[i]] *
           cos(theta_i - PHI_STIM[which_pop[i]] * M_PI / 180.0));
    }
  }

  if (step == (int) (T_STIM[1] / DT)) {
    if (VERBOSE)
      std::cout << " STIM OFF" << std::endl;
    for (int i = 0; i < N; i++)
      ff_inputs[i] = Iext_scaled[which_pop[i]];
  }

  if (IF_FF_NOISE) {
    for (int i = 0; i < N; i++)
      ff_inputs[i] += sqrt(VAR_FF[which_pop[i]] / 1000.0) * white(gen);
  }
}

void updateVolts(){
  int pres_pop=0;

  for (int i = 0; i < N; ++i) {
    pres_pop = which_pop[i];
    volts[i] *= EXP_DT_TAU_MEM[pres_pop];
    // volts[i] += DT_TAU_MEM[pres_pop] * (net_inputs[i] + V_LEAK);
    volts[i] += DT * (net_inputs[i] + V_LEAK); // This gets me the correct mf rates
  }

  // float RK1=0, RK2=0;
  // for (int i = 0; i < N; ++i) {
  //   pres_pop = which_pop[i];
  //   RK1 = -(volts[i] - V_LEAK) / TAU_MEM[pres_pop] + net_inputs[i] ;
  //   RK2 = -(volts[i] - V_LEAK + DT * RK1) / TAU_MEM[pres_pop] + net_inputs_RK2[i] ;
  //   volts[i] = volts[i] + DT / 2.0 * ( RK1 + RK2 ) ;
  // }

}

void updateRecInputs(){
  int pres_pop=0, post_pop=0;

  for (int j = 0; j < N; ++j) // presynaptic
    if (spikes[j] == 1) {
      pres_pop = which_pop[j];

      for (size_t i = colptr[j]; i < colptr[j + 1]; ++i) { // postsynaptic
        post_pop = which_pop[indices[i]];
        if (IF_STP && pres_pop == 0 && post_pop == 0)
          inputs[pres_pop][indices[i]] += A_stp[j] * Jab_scaled[pres_pop + N_POP * post_pop];
        else
          inputs[pres_pop][indices[i]] += Jab_scaled[pres_pop + N_POP * post_pop];
      }
    }

  if (IF_NMDA) {
    for (int j = 0; j < Na[0]; ++j) // presynaptic
      if (spikes[j] == 1) {
        for (size_t i = colptr[j]; i < colptr[j + 1]; ++i) { // postsynaptic
          post_pop = which_pop[indices[i]];
          if (IF_STP && post_pop == 0)
            inputs_NMDA[indices[i]] += A_stp[j] * Jab_NMDA[post_pop];
          else
            inputs_NMDA[indices[i]] += Jab_NMDA[post_pop];
        }
      }
  }
}

void updateNetInputs(){

  for(int i = 0; i < N; ++i)
    net_inputs[i] = ff_inputs[i];

  for(int i = 0; i < N_POP; ++i) {
    for (int j = 0; j < N; ++j) {
      net_inputs[j] += inputs[i][j];
      inputs[i][j] *= EXP_DT_TAU_SYN[i];
    }
  }
  if (IF_NMDA) {
    for (int j = 0; j < N; ++j) {
      net_inputs[j] += inputs_NMDA[j];
      inputs_NMDA[j] *= EXP_DT_TAU_NMDA[which_pop[j]];
    }
  }
}

void updateSpikes(int step){
  for (int i = 0; i < N; ++i)
    if (volts[i] >= V_THRESH) {
      volts[i] = V_REST;
      spikes[i] = 1.0;

      // if (IF_RK2) {
      //   dt2 = DT * (V_THRESH - vold[i]) / (volts[i] - vold[i] );
      //   volts[i] = (volts[i] - V_THRESH)
      //     * (1.0 + DT_TAU_MEM[pop] * (vold[i] - V_REST)
      //        / (volts[i] - vold[i])) + V_REST;
      // }

      if (IF_STP && i < Na[0]) updateStp(i, step);

    } else {
      spikes[i] = 0.0;
    }

  for (int i = 0; i < N; ++i)
    rates[i] += spikes[i];
}

void updateStp(int i, int step){
  // This is the Mato & Hansel stp model
  float ISI = step * DT - spike_times[i];
  spike_times[i] = step * DT;

  int pre_pop = which_pop[i];

  A_stp[i] = u_stp[i] * x_stp[i];
  u_stp[i] = u_stp[i] * exp(-ISI / TAU_FAC[pre_pop]) + USE[pre_pop] * (1.0 - u_stp[i] * exp(-ISI / TAU_FAC[pre_pop]));
  x_stp[i] = x_stp[i] * (1.0 - u_stp[i]) * exp(-ISI / TAU_REC[pre_pop]) + 1.0 - exp(-ISI / TAU_REC[pre_pop]);

  // // This is the Mongillo et al. 2012 model
  // if (dist(gen) < U) { // calcium binding with proba U
  //   u_stp[i] = 1.0;
  //   if (x_stp[i] == 1.0) { // neurotransmitter release if available
  //     A_stp[i] = 1.0; // spike
  //     x_stp[i] = 0.0;
  //   }
  // }

  // // In between spikes
  // if (u_stp[i]==1.0 && dist(gen) < 1.0 / TAU_FAC[pre_pop])  // calcium unbinds with proba 1/TAU_FAC
  //   u_stp[i] = 0.0;
  // if (x_stp[i]==0.0 && dist(gen) < 1.0 / TAU_FAC[pre_pop])  // neurotransmitter refill with proba 1/TAU_REC
  //   x_stp[i] = 1.0;
}

void printParam(){
  std::cout << "N_POP " << N_POP;

  std::cout << " N " << N << " Na ";
  for(int i=0; i < N_POP; i++)
    std::cout << Na[i] << " ";

  std::cout << " cNa ";
  for(int i=0; i < N_POP+1; i++)
    std::cout << cNa[i] << " ";
  std::cout << std::endl;

  std::cout << "K " << K << " Ka ";
  for(int i=0; i < N_POP; i++)
    std::cout << Ka[i] << " ";
  std::cout << std::endl;

  std::cout << "Jab ";
  for(int i=0; i < N_POP; i++)
    for(int j=0; j < N_POP; j++)
      std::cout << Jab[j + i * N_POP] << " ";
  // std::cout << std::endl;

  std::cout << "Iext ";
  for(int i=0; i < N_POP; i++)
    std::cout << Iext[i] << " ";
  std::cout << std::endl;

  std::cout << "Iext_scaled ";
  for (int i = 0; i < N_POP; i++) {
    std::cout << Iext_scaled[i] << " ";
  }
  // std::cout << std::endl;

  std::cout << "Jab_scaled ";
  for(int i=0; i < N_POP; i++)
    for(int j=0; j < N_POP; j++)
      std::cout << Jab_scaled[j + i * N_POP] << " ";
  std::cout << std::endl;
}

void runSimul(){

  const float dum = 1000.0 / T_WINDOW;

  initNetwork();

  if(IF_LOAD_MAT)
    getSparseMatCSC(colptr, indices);
  else
    genSparseMatCSC(colptr, indices);

  if (VERBOSE)
    printParam();

  std::ofstream (DATA_PATH + "/rates.txt", std::ios::trunc).close();
  std::ofstream (DATA_PATH + "/spikes.txt", std::ios::trunc).close();
  std::ofstream (DATA_PATH + "/inputsE.txt", std::ios::trunc).close();
  std::ofstream (DATA_PATH + "/inputsI.txt", std::ios::trunc).close();
  std::ofstream (DATA_PATH + "/volts.txt", std::ios::trunc).close();

  std::ofstream (DATA_PATH + "/x_stp.txt", std::ios::trunc).close();
  std::ofstream (DATA_PATH + "/u_stp.txt", std::ios::trunc).close();
  std::ofstream (DATA_PATH + "/A_stp.txt", std::ios::trunc).close();

  std::ofstream ratesFile(DATA_PATH + "/rates.txt", std::ios::app | std::ios::binary);
  std::ofstream spikesFile(DATA_PATH + "/spikes.txt", std::ios::app | std::ios::binary);
  std::ofstream inputsEfile(DATA_PATH + "/inputsE.txt", std::ios::app | std::ios::binary);
  std::ofstream inputsIfile(DATA_PATH + "/inputsI.txt", std::ios::app | std::ios::binary);
  std::ofstream voltsFile(DATA_PATH + "/volts.txt", std::ios::app | std::ios::binary);

  std::ofstream xstpFile(DATA_PATH + "/x_stp.txt", std::ios::app | std::ios::binary);
  std::ofstream ustpFile(DATA_PATH + "/u_stp.txt", std::ios::app | std::ios::binary);
  std::ofstream AstpFile(DATA_PATH + "/A_stp.txt", std::ios::app | std::ios::binary);

  int N_STEPS = (int) (DURATION/DT);
  int N_STEADY = (int) T_STEADY / DT;
  int N_WINDOW = (int) T_WINDOW / DT;

  std::cout << "Running Simulation" << std::endl;
  for(int step = 0; step < N_STEPS; step += 1) {

    updateVolts();
    updateSpikes(step); // must come before updateRecInputs in this implementation
    updateFFinputs(step); // must come before updateNetInputs in this implementation
    updateRecInputs(); // must come before updateNetInputs in this implementation
    updateNetInputs();

    if(step % N_WINDOW == 0 && step >= N_STEADY) {
      for(int i=0; i<N; ++i)
        rates[i] *= dum;

      if (VERBOSE) {
        std::cout << std::setprecision(3);
        std::cout << "time " << step * DT << "s";

        std::cout << "| Rates ";
        for (int i = 0; i < N_POP; ++i)
          std::cout << popMean(rates, cNa[i], cNa[i + 1]) << " Hz ";

        std::cout << "| Spike count ";
        for (int i = 0; i < N_POP; ++i)
          std::cout << popMean(spikes, cNa[i], cNa[i + 1]) * 1000.0 / DT << " ";
        std::cout << std::flush;
        std::cout << "\r";
      }

      saveArrayToFile(spikesFile, spikes, N);
      saveArrayToFile(ratesFile, rates, N);

      saveArrayToFile(voltsFile, volts, N);

      saveArrayToFile(inputsEfile, inputs[0], N);
      saveArrayToFile(inputsIfile, inputs[1], N);

      if (IF_STP) {
        saveArrayToFile(xstpFile, x_stp, Na[0]);
        saveArrayToFile(ustpFile, u_stp, Na[0]);
        saveArrayToFile(AstpFile, A_stp, Na[0]);
      }

      for(int i=0; i<N; ++i)
        rates[i] = 0.0 ;
    } // end window
  } // end for

  spikesFile.close();
  ratesFile.close();
  inputsEfile.close();
  inputsIfile.close();
  voltsFile.close();

  xstpFile.close();
  ustpFile.close();
  AstpFile.close();

  std::cout << "Done" << std::endl;
}
