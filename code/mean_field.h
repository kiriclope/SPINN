#ifndef __MFRATES__
#define __MFRATES__

#include <iostream>
// #include <eigen3/Eigen/Dense>
/* #include <eigen3/unsupported/Eigen/NonLinearOptimization> */
#include "../../libs/eigen/eigen-3.4.0/unsupported/Eigen/NonLinearOptimization"

using namespace Eigen ;

void mean_field_rates() { 

  bool balance_condition = false ; 
  float JeJi, IeIi ;
  
  if(n_pop==2)
    if(IF_STP==0 ) {
      IeIi = ( ext_inputs[0] / J[0] ) / ( ext_inputs[1] / J[2] ) ; 
      JeJi = abs( J[1] / J[0] ) / abs( J[3] / J[2] ) ; 
      balance_condition = IeIi>JeJi && JeJi>1.0 ; 
    }
    else{
      IeIi = abs( ext_inputs[0] / J[1] ) / abs( ext_inputs[1] / J[3] ) ; 
      balance_condition = IeIi>1.0 ; 
    } 
  else
    balance_condition = 1.0 ;
  
  cout << "###############################################" << endl ;
  
  if( balance_condition) {
    
    if(IF_STP==0)
      cout << "BALANCE CONDITIONS VALID: Ie/Ii > Je/Ji > 1 " << endl ; 
    else 
      cout << "BALANCE CONDITIONS VALID: Ie/Jei > Ii/Jii  " << endl ; 
    
    mf_rates = new float [n_pop]() ;

    MatrixXd A(n_pop,n_pop) ;
    VectorXd b(n_pop) ;
    
    for(int i=0;i<n_pop;i++) {
      /* b(i) = GAIN*ext_inputs[i]*m0*1000.  - Vth / TAU_MEM[i] / sqrt_Ka[0] ;  */
      b(i) = GAIN*ext_inputs[i]*m0*1000. ;
      for(int j=0;j<n_pop;j++) 
	A(i,j) = GAIN*J[j+i*n_pop] ; 
    }
    
    VectorXd x(n_pop) ; 
    x = A.colPivHouseholderQr().solve(-b) ; 
    
    for(int i=0;i<n_pop;i++)
      mf_rates[i] = x(i) ;
  }
  else {
    cout << "ERROR BALANCE CONDITIONS INVALID: " ;
    if(IF_STP==0) {
      if (JeJi<1)
	cout << "Je/Ji < 1 " << endl ;
      else
	cout << "Ie/Ii < Je/Ji " << endl ; 
    }
    else
      cout << "Ie/Jei < Ii/Jii  " << endl ; 
  }  
  
}

#endif
