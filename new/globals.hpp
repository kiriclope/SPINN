#include <cmath>

#define N (int) 1000
#define N_POP (int) 1
#define K (float) 1.0

#define I0 (float) 10.0
#define J0 (float) -2.0
#define J1 (float) 0.4

#define DT (float) 0.0005
#define DURATION (float) 10.0
#define T_WINDOW (float) 0.05

#define THRESH (float) 1.0

#define TAU_SYN (float) 0.002
#define DT_TAU_SYN (float) DT / TAU_SYN
#define EXP_DT_TAU_SYN (float) exp(-DT / TAU_SYN)

#define TAU_RATE (float) 0.020
#define DT_TAU_RATE (float) DT / TAU_RATE
#define EXP_DT_TAU_RATE (float) exp(-DT / TAU_RATE)
