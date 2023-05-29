#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>

#include "globals.hpp"
#include "utils.hpp"
#include "sparse_mat.hpp"
#include "lif_network.hpp"

void LifNetwork::initNetwork() {

  // volts = generateGaussianVector<float>(N, 2.0, 0.5);
  for(int i=0; i < N_POP; i++)
    for(int j=0; j < N_POP; j++)
      Jab_scaled[j + i * N_POP] = GAIN * Jab[j + i * N_POP] * (V_THRESH - V_REST) / TAU_SYN[j] / sqrt(K);

  // Here (V_THRESH - V_REST) is important to get rates from the balanced linear eqs m = I/J;

  for(int i=0; i < N_POP; i++)
    Iext_scaled[i] = GAIN * Iext[i] * sqrt(K) * (V_THRESH - V_REST);

  for(int i=0; i<N; i++)
    ff_inputs[i] = Iext_scaled[which_pop[i]];
}

void LifNetwork::updateFFinputs(int step) {

  float theta_i = 0;

  if (step == (int) (T_STIM[0] / DT)) {
    std::cout << " STIM ON" << std::endl;
    for (int i = 0; i < N; i++) {
      theta_i = (2.0 * M_PI * i) / (float) Na[which_pop[i]];

      ff_inputs[i] = Iext_scaled[which_pop[i]]
        + A_STIM[which_pop[i]] * sqrt(K)
        * (1.0 + KAPPA_STIM[which_pop[i]] *
           cos(theta_i - PHI_STIM[which_pop[i]] * M_PI / 180.0));
    }
  }

  if (step == (int) (T_STIM[1] / DT)) {
    std::cout << " STIM OFF" << std::endl;
    for (int i = 0; i < N; i++)
      ff_inputs[i] = Iext_scaled[which_pop[i]];
  }

  if (IF_FF_NOISE) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<float> dist(0.0, 1.0);
    for (int i = 0; i < N; i++)
      ff_inputs[i] += dist(gen) * sqrt(K);
  }
}

void LifNetwork::updateVolts(){
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

void LifNetwork::updateInputs(){
  int pres_pop=0, post_pop=0;

  for (int j = 0; j < N; ++j) // presynaptic
    if (spikes[j] == 1) {
      pres_pop = which_pop[j];

      for(unsigned long i = colptr[j]; i < colptr[j + 1]; ++i) { // postsynaptic
        post_pop = which_pop[indices[i]];
        if (IF_STP && pres_pop==0 && post_pop==0)
          inputs[pres_pop][indices[i]] += A_stp[j] * Jab_scaled[pres_pop + N_POP * post_pop];
        else
          inputs[pres_pop][indices[i]] += Jab_scaled[pres_pop + N_POP * post_pop];
      }
    }

  for(int i = 0; i < N; ++i)
    net_inputs[i] = ff_inputs[i];

  for(int i = 0; i < N_POP; ++i) {
    for (int j = 0; j < N; ++j) {
      net_inputs[j] += inputs[i][j];
      inputs[i][j] *= EXP_DT_TAU_SYN[i];
    }
  }
}

void LifNetwork::updateSpikes(int step){
  for (int i = 0; i < N; ++i)
    if (volts[i] >= V_THRESH) {
      volts[i] = V_REST;
      spikes[i] = 1.0;

      if(IF_STP)
        updateStp(i, step);

    } else {
      spikes[i] = 0.0;
    }

  for (int i = 0; i < N; ++i)
    rates[i] += spikes[i];
}

void LifNetwork::updateStp(int i, int step){
  // This is the Mato & Hansel stp model
  float ISI = step * DT - spike_times[i];
  spike_times[i] = step * DT;

  int pre_pop = which_pop[i];

  u_stp[i] =
      u_stp[i] * exp(- ISI / TAU_FAC[pre_pop]) +
      USE[pre_pop] * (1.0 - u_stp[i] * exp(-ISI / TAU_FAC[pre_pop]));

  x_stp[i] =
      x_stp[i] * (1.0 - u_stp[i]) * exp(-ISI / TAU_REC[pre_pop]) +
      1.0 - exp(- ISI / TAU_REC[pre_pop]);

  A_stp[i] = u_stp[i] * x_stp[i];
}

void LifNetwork::printParam(){
  std::cout << "LIF NETWORK ";
  std::cout << "N_POP " << N_POP << std::endl;

  std::cout << "N " << N << " Na ";
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
  std::cout << std::endl;

  std::cout << "Iext ";
  for(int i=0; i < N_POP; i++)
    std::cout << Iext[i] << " ";
  std::cout << std::endl;

  std::cout << "Iext_scaled ";
  for (int i = 0; i < N_POP; i++) {
    std::cout << Iext_scaled[i] << " ";
  }
  std::cout << std::endl;

  std::cout << "Jab_scaled ";
  for(int i=0; i < N_POP; i++)
    for(int j=0; j < N_POP; j++)
      std::cout << Jab_scaled[j + i * N_POP] << " ";
  std::cout << std::endl;
}

void LifNetwork::runSimul(){

  const float dum = 1.0 / T_WINDOW;

  std::cout << "Initializing Network";
  initNetwork();
  std::cout << " Done" << std::endl;

  std::cout << "Generating Sparse Matrix";
  genSparseMatCSC(colptr, indices);
  std::cout << " Done" << std::endl;

  printParam();

  std::ofstream ratesFile("rates.txt");
  std::ofstream spikesFile("spikes.txt");
  std::ofstream inputsEfile("inputsE.txt");
  std::ofstream inputsIfile("inputsI.txt");
  std::ofstream voltsFile("volts.txt");

  std::cout << "Running Simulation" << std::endl;
  for(int step = 0; step < DURATION/DT; step += 1) {

    updateVolts();
    updateSpikes(step); // must come before updateInputs in this implementation
    updateFFinputs(step);
    updateInputs();

    if(step % (int) (T_WINDOW / DT) == 0) {
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
          std::cout << popMean(spikes, cNa[i], cNa[i + 1]) / DT << " ";
        std::cout << std::flush;
        std::cout << "\r";

        saveArrayToFile(spikesFile, spikes, N);
        saveArrayToFile(ratesFile, rates, N);

        saveArrayToFile(voltsFile, volts, N);

        saveArrayToFile(inputsEfile, inputs[0], N);
        saveArrayToFile(inputsIfile, inputs[1], N);

        for(int i=0; i<N; ++i)
          rates[i] = 0.0 ;
      }

    }
  }

  spikesFile.close();
  ratesFile.close();
  inputsEfile.close();
  inputsIfile.close();
  voltsFile.close();

  std::cout << "Done" << std::endl;
}
