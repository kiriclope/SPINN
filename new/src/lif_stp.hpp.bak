#include "globals.hpp"
#include "lif_stp.hpp"

void matoModel(int i, float ISI) {
  pre_pop = which_pop[i];

  u_stp[i] =
      u_stp[i] * exp(-ISI / TAU_FAC[pre_pop]) +
      USE[pre_pop] * (1. - u_stp[i] * exp(-ISI / TAU_FAC[pre_pop]));

  x_stp[i] =
      x_stp[i] * (1. - u_stp[i]) * exp(-ISI / TAU_REC[pre_pop]) +
      1. - exp(-ISI / TAU_REC[pre_pop]);
  A_stp[i] = u_stp[i] * x_stp[i];
}
