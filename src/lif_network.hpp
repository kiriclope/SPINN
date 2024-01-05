#ifndef NETWORK_HPP
#define NETWORK_HPP

#include "globals.hpp"

extern float *rates;
extern float *volts;
extern float *spikes;
extern float *spike_times;

extern float *ff_inputs;
extern float **inputs;
extern float *inputs_NMDA;
extern float *net_inputs;

extern float *Jab_scaled;
extern float *Jab_NMDA;

extern float *Iext_scaled;

extern unsigned long *colptr;
extern int *indices;

extern float **x_stp;
extern float **u_stp;
extern float **A_stp;

void init_lif();
void free_lif();

void printParam();
void initNetwork();
void updateVolts();
void updateFFinputs(int step);
void updateRecInputs();
void updateNetInputs();
void updateSpikes(int step);
void updateStp(int i, int post_pop, float isi);
void runSimul();

void odr_stimuli(float* &ff_inputs, int FLAG);
void dual_task_stimuli(float* &ff_inputs, int FLAG);

#endif
