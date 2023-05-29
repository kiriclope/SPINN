#ifndef GLOBALS_HPP
#define GLOBALS_HPP
#include <vector>
#include <string>

// Declare all variables as extern
extern int VERBOSE;
extern int N;
extern int N_POP;
extern float K;

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

extern std::vector<float> TAU_SYN;
extern float* DT_TAU_SYN;
extern float* EXP_DT_TAU_SYN;

extern float V_THRESH;
extern float V_REST;
extern float V_LEAK;

extern std::vector<float> TAU_MEM;
extern float* DT_TAU_MEM;
extern float* EXP_DT_TAU_MEM;

extern std::string PROBA;
extern std::vector<float> KAPPA;

extern int IF_STP;
extern std::vector<float> TAU_REC;
extern std::vector<float> TAU_FAC;
extern std::vector<float> USE;

extern std::vector<float> T_STIM ;
extern std::vector<float> A_STIM ;
extern std::vector<float> PHI_STIM ;
extern std::vector<float> KAPPA_STIM ;

extern int IF_FF_NOISE;
extern std::vector<float> VAR_FF;

void loadConfig(std::string configname);
#endif
