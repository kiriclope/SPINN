#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <yaml-cpp/yaml.h>
#include <random>
// #define ARMA_DONT_PRINT_FAST_MATH_WARNING
// #include <armadillo>

#include "globals.hpp"

std::normal_distribution<float> white(0.0, 1.0);
std::uniform_real_distribution<float> unif(0.0, 1.0);

std::string DATA_PATH;
std::string MAT_PATH;

int CHECK_BISTABILITY;
std::vector<int> BUMP_SWITCH;

int VERBOSE;
int N;
int N_POP;
float K;

std::vector<float> FRAC;
int* Na ;
int* cNa ;
float* Ka ;
int* which_pop ;

float GAIN;
std::vector<float> Jab ;
std::vector<float> Iext ;

float DT;
float DURATION;
float T_WINDOW;
float T_STEADY;

std::vector<float> TAU_SYN;
float* DT_TAU_SYN;
float* EXP_DT_TAU_SYN;

float V_THRESH;
float V_REST;
float V_LEAK;

std::vector<float> TAU_MEM;
float* DT_TAU_MEM;
float* EXP_DT_TAU_MEM;

int IF_LOAD_MAT;
int IF_SAVE_MAT;
int IF_SAVE_DATA;
float T_SAVE;

std::vector<std::string> PROBA;
// std::string PROBA;
std::vector<float> KAPPA;

// LOW RANK
// arma::mat ksi;
int LR_RANK;
int LR_SEED;
int LR_LOAD;

std::vector<float> LR_MEAN;
std::vector<float> LR_STD;
std::vector<float> LR_RHO;
std::vector<float> LR_FF_RHO;

std::vector<float> ksi_0(N);
std::vector<float> ksi_1(N);
std::vector<float> ksi_2(N);

int IF_NMDA;
std::vector<float> TAU_NMDA;
std::vector<float> R_NMDA;
float* EXP_DT_TAU_NMDA;

int IF_STP;
std::vector<float> TAU_REC;
std::vector<float> TAU_FAC;
std::vector<float> USE;

std::vector<float> T_STIM ;
int* N_STIM ;
std::vector<float> A_STIM ;
std::vector<float> STD_STIM ;
std::vector<float> PHI_STIM ;
std::vector<float> KAPPA_STIM ;

std::vector<float> T_DIST ;
int* N_DIST ;
std::vector<float> A_DIST ;
std::vector<float> STD_DIST ;
std::vector<float> PHI_DIST ;
std::vector<float> KAPPA_DIST ;

int IF_FF_NOISE;
std::vector<float> STD_FF;

int IF_FF_CORR;
std::vector<float> A_CORR;
std::vector<float> CORR_FF;
double phi0;

void loadConfig(std::string configname){

  YAML::Node config = YAML::LoadFile(configname);

  DATA_PATH = config["DATA_PATH"].as<std::string>();
  MAT_PATH = config["MAT_PATH"].as<std::string>();

  IF_LOAD_MAT = config["IF_LOAD_MAT"].as<int>();
  IF_SAVE_MAT = config["IF_SAVE_MAT"].as<int>();
  
  IF_SAVE_DATA = config["IF_SAVE_DATA"].as<int>();
  T_SAVE = config["T_SAVE"].as<float>();
  
  VERBOSE = config["verbose"].as<int>();
  CHECK_BISTABILITY = config["CHECK_BISTABILITY"].as<int>();
  BUMP_SWITCH = config["BUMP_SWITCH"].as<std::vector<int>>();
  
  // Assign variables from configuration file
  N = config["N"].as<int>();
  N_POP = config["N_POP"].as<int>();
  K = config["K"].as<float>();

  FRAC = config["FRAC"].as<std::vector<float>>();

  Na = new int[N_POP]() ;
  Ka = new float[N_POP]() ;
  cNa = new int[N_POP+1]() ;
  which_pop = new int[N]() ;

  for(int i=0; i<N_POP; ++i) {
    Na[i] = (int) (FRAC[i] * N);
    // Ka[i] = K;
    Ka[i] = FRAC[i] * K;
  }

  cNa[0] = 0;
  cNa[1] = Na[0] ;
  cNa[2] = Na[0] + Na[1];

  for(int i = 0; i < N; i++) {
    for (int i_pop = 0; i_pop < N_POP; ++i_pop)
      if (i >= cNa[i_pop] && i < cNa[i_pop + 1])
        which_pop[i] = i_pop;
  }

  GAIN = config["GAIN"].as<float>();
  Jab = config["Jab"].as<std::vector<float>>();
  Iext = config["Iext"].as<std::vector<float>>();

  DT = config["DT"].as<float>();
  DURATION = config["DURATION"].as<float>();
  T_WINDOW = config["T_WINDOW"].as<float>();
  T_STEADY = config["T_STEADY"].as<float>();

  TAU_SYN = config["TAU_SYN"].as<std::vector<float>>();

  DT_TAU_SYN = new float[N_POP]();
  EXP_DT_TAU_SYN = new float[N_POP]();

  for(int i=0; i<N_POP; ++i) {
    DT_TAU_SYN[i] = DT / TAU_SYN[i];
    EXP_DT_TAU_SYN[i] = std::exp(-DT / TAU_SYN[i]);
  }

  V_THRESH = config["V_THRESH"].as<float>();
  V_REST = config["V_REST"].as<float>();
  V_LEAK = config["V_LEAK"].as<float>();

  TAU_MEM = config["TAU_MEM"].as<std::vector<float>>();

  DT_TAU_MEM = new float[N_POP]();
  EXP_DT_TAU_MEM = new float[N_POP]();

  for(int i=0; i<N_POP; ++i) {
    DT_TAU_MEM[i] = DT / TAU_MEM[i];
    EXP_DT_TAU_MEM[i] = std::exp(-DT / TAU_MEM[i]);
  }

  PROBA = config["PROBA"].as<std::vector<std::string>>();
  KAPPA = config["KAPPA"].as<std::vector<float>>();

  IF_NMDA = config["IF_NMDA"].as<int>();
  TAU_NMDA = config["TAU_NMDA"].as<std::vector<float>>();
  R_NMDA = config["R_NMDA"].as<std::vector<float>>();

  EXP_DT_TAU_NMDA = new float[N_POP]();
  for(int i=0; i<N_POP; ++i)
    EXP_DT_TAU_NMDA[i] = std::exp(-DT / TAU_NMDA[i]);

  IF_STP = config["IF_STP"].as<int>();
  TAU_REC = config["TAU_REC"].as<std::vector<float>>();
  TAU_FAC = config["TAU_FAC"].as<std::vector<float>>();
  USE = config["USE"].as<std::vector<float>>();

  T_STIM = config["T_STIM"].as<std::vector<float>>();
  N_STIM = new int[N_POP]();
  for(int i=0; i<N_POP; ++i)
    N_STIM[i] = (int) ((T_STIM[i] + T_STEADY) / DT);
  
  A_STIM = config["A_STIM"].as<std::vector<float>>();
  STD_STIM = config["STD_STIM"].as<std::vector<float>>();
  PHI_STIM = config["PHI_STIM"].as<std::vector<float>>();
  KAPPA_STIM = config["KAPPA_STIM"].as<std::vector<float>>();

  T_DIST = config["T_DIST"].as<std::vector<float>>();
  N_DIST = new int[N_POP]();
  for(int i=0; i<N_POP; ++i)
    N_DIST[i] = (int) ((T_DIST[i] + T_STEADY) / DT);
  
  A_DIST = config["A_DIST"].as<std::vector<float>>();
  STD_DIST = config["STD_DIST"].as<std::vector<float>>();
  PHI_DIST = config["PHI_DIST"].as<std::vector<float>>();
  KAPPA_DIST = config["KAPPA_DIST"].as<std::vector<float>>();
  
  IF_FF_NOISE = config["IF_FF_NOISE"].as<int>();
  STD_FF = config["STD_FF"].as<std::vector<float>>();

  IF_FF_CORR = config["IF_FF_CORR"].as<int>();
  CORR_FF = config["CORR_FF"].as<std::vector<float>>();
  A_CORR = config["A_CORR"].as<std::vector<float>>();

  
  LR_RANK = config["LR_RANK"].as<int>();
  LR_SEED = config["LR_SEED"].as<int>();
  LR_LOAD = config["LR_LOAD"].as<int>();
  
  LR_MEAN = config["LR_MEAN"].as<std::vector<float>>();
  LR_STD = config["LR_STD"].as<std::vector<float>>();
  LR_RHO = config["LR_RHO"].as<std::vector<float>>();
  LR_FF_RHO = config["LR_FF_RHO"].as<std::vector<float>>();

  if(PHI_STIM[0] == 180.0)
    LR_FF_RHO[0] *= -1;
  
}

void ensureDirExists(std::string &path) {
    if (!std::filesystem::exists(path)) {
      std::filesystem::create_directories(path);
    }
}
