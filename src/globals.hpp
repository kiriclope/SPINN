#ifndef GLOBALS_HPP
#define GLOBALS_HPP
#include <vector>
#include <string>
#include <filesystem>
#include <string>
#include <random>
//#define ARMA_DONT_PRINT_FAST_MATH_WARNING
//#include <armadillo>

extern std::normal_distribution<float> white;
extern std::uniform_real_distribution<float> unif;

// Declare all variables as extern
extern std::string DATA_PATH;
extern std::string MAT_PATH;

extern int CHECK_BISTABILITY;
extern std::vector<int> BUMP_SWITCH;

extern int VERBOSE;
extern int N;
extern int N_POP;
extern float K;

extern float M0;
extern std::vector<float> FRAC;
extern int* Na ;
extern int* cNa ;
extern float* Ka ;
extern int* which_pop ;

extern float GAIN;
extern std::vector<float> Iext ;
extern std::vector<float> Jab ;

extern float DT;
extern float DURATION;
extern float T_WINDOW;
extern float T_STEADY;

extern std::vector<float> TAU_SYN;
extern float* DT_TAU_SYN;
extern float* EXP_DT_TAU_SYN;

extern float V_THRESH;
extern float V_REST;
extern float V_LEAK;
extern std::vector<float> V_REV;

extern int IF_COND_BASE;
extern int IF_THRESH_DYN;
extern float DELTA_THRESH;
extern std::vector<float> TAU_AREF;
extern std::vector<float> TAU_REF;
extern float* DT_TAU_REF;
extern float* EXP_DT_TAU_REF;

extern std::vector<float> TAU_MEM;
extern float* DT_TAU_MEM;
extern float* EXP_DT_TAU_MEM;

extern int IF_SAVE_DATA;
extern float T_SAVE;

extern int IF_LOAD_MAT;
extern int IF_SAVE_MAT;
extern std::vector<std::string> PROBA;
// extern std::string PROBA;
extern std::vector<float> KAPPA;

// extern arma::mat ksi;
extern int LR_RANK;
extern int LR_SEED;
extern int LR_LOAD;

extern std::vector<float> LR_MEAN;
extern std::vector<float> LR_STD;
extern std::vector<float> LR_RHO;

extern std::vector<float> ksi_0;
extern std::vector<float> ksi_1;
extern std::vector<float> ksi_2;

extern int IF_NMDA;
extern std::vector<float> TAU_NMDA;
extern std::vector<float> R_NMDA;
extern float* EXP_DT_TAU_NMDA;

extern int IF_STP;
extern std::vector<float> TAU_REC;
extern std::vector<float> TAU_FAC;
extern std::vector<float> USE;

extern std::vector<float> T_STIM ;
extern int* N_STIM ;
extern std::vector<float> A_STIM ;
extern std::vector<float> STD_STIM ;
extern std::vector<float> PHI_STIM ;
extern std::vector<float> KAPPA_STIM ;

extern std::vector<float> T_DIST ;
extern int* N_DIST ;
extern std::vector<float> A_DIST ;
extern std::vector<float> STD_DIST ;
extern std::vector<float> PHI_DIST ;
extern std::vector<float> KAPPA_DIST ;

extern int IF_FF_NOISE;
extern std::vector<float> STD_FF;

extern int IF_FF_CORR;
extern std::vector<float> CORR_FF;
extern std::vector<float> A_CORR;
extern double phi0;

void loadConfig(std::string configname);
void ensureDirExists(std::string &path);

#endif
