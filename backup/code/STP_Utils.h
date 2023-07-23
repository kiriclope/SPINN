#ifndef __STPUTILS__
#define __STPUTILS__

//STP Parameters 
#define IF_STP 1
const int STP[16] = {1,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,0} ;

const double Trec[4] = {0.2E3,  0.2E3, 0.2E3, 0.3E3} ; 
const double Tfac[4] = {3.0E3,  0.6E3, 0.6E3, 1.E3} ; 
const double Use[4] = { 0.07,    0.05,   0.05,  0.3} ; 

void Markram98(double &u, double& x, double ISI, int whichPop) {
  u = u * exp(- ISI / Tfac[whichPop] ) + Use[whichPop] * ( 1. - u * exp(- ISI / Tfac[whichPop] ) ) ; 
  x = x * ( 1. - u ) * exp(- ISI / Trec[whichPop] ) + 1. - exp(- ISI / Trec[whichPop] ) ; 
}

void Mongillo08(double &u, double &x, int whichPop) {
  x += DT * ( 1.0 - x ) / Trec[whichPop] ; 
  u += DT * ( Use[whichPop] - u ) / Tfac[whichPop] ; 
}

void STP_Rate(double &u, double &x, double &r, int whichPop) {
  x += DT * ( ( 1.0 - x ) / Trec[whichPop]  + Use[whichPop] * (1.0 - u) * r ) ; 
  u += DT * ( ( Use[whichPop] - u ) / Tfac[whichPop]  - u * x * r ) ; 
}

#define IF_Mongillo12 1
#define IF_MATO 0

#endif
