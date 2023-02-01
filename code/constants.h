#ifndef __CONSTANTS__ 
#define __CONSTANTS__ 

//////////////////////////////////////////
// Simulation globals 
//////////////////////////////////////////
#define DT .01 
#define DURATION 2.00E3 // 10.E3 // 

#define TIME_STEADY 1.0E3 //  2.0E3 // 10.E3 // 
#define TIME_WINDOW 0.0125E3 // 1.00E3 // 10.E3 // 
#define TIME_REC 1.0E3 // 1.00E3 // 10.E3 // 

#define IF_EULER 1 
#define IF_RK2 0 

////////////////////////////////////////// 
// Network globals 
////////////////////////////////////////// 
const int eps[2] = {1,-1} ; 
const double Trate[2] = {20.0, 10.0} ; 

#define J0 -1.0 // for one population 
#define I0 0.5 // for one population 
#define Tsyn0 1.0 // for one population 

#define GAIN 1.0 
#define m0 0.1 

////////////////////////////////////////// 
// Loop globals 
////////////////////////////////////////// 
#define IF_TRIALS 0 
#define TRIAL_ID 6

#endif
