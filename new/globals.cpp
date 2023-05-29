#include <vector>
#include <string>
#include <cmath>
#include <yaml-cpp/yaml.h>

#include "globals.hpp"

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

std::vector<float> TAU_SYN;
float* DT_TAU_SYN;
float* EXP_DT_TAU_SYN;

float V_THRESH;
float V_REST;
float V_LEAK;

std::vector<float> TAU_MEM;
float* DT_TAU_MEM;
float* EXP_DT_TAU_MEM;

std::string PROBA;
std::vector<float> KAPPA;

int IF_STP;
std::vector<float> TAU_REC;
std::vector<float> TAU_FAC;
std::vector<float> USE;

std::vector<float> T_STIM ;
std::vector<float> A_STIM ;
std::vector<float> PHI_STIM ;
std::vector<float> KAPPA_STIM ;

int IF_FF_NOISE;
std::vector<float> VAR_FF;

void loadConfig(std::string configname){
  YAML::Node config = YAML::LoadFile(configname);

  VERBOSE = config["verbose"].as<int>();

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

  PROBA = config["PROBA"].as<std::string>();
  KAPPA = config["KAPPA"].as<std::vector<float>>();

  IF_STP = config["IF_STP"].as<int>();
  TAU_REC = config["TAU_REC"].as<std::vector<float>>();
  TAU_FAC = config["TAU_FAC"].as<std::vector<float>>();
  USE = config["USE"].as<std::vector<float>>();

  T_STIM = config["T_STIM"].as<std::vector<float>>();
  A_STIM = config["A_STIM"].as<std::vector<float>>();
  PHI_STIM = config["PHI_STIM"].as<std::vector<float>>();
  KAPPA_STIM = config["KAPPA_STIM"].as<std::vector<float>>();

  IF_FF_NOISE = config["IF_FF_NOISE"].as<int>();
  VAR_FF = config["VAR_FF"].as<std::vector<float>>();

}
