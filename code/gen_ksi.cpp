#include <iostream> 
#include <iomanip> 
#include <sstream>
#include <fstream>
#include <cmath> 
#include <random>
#include <ctime> 
#include <sys/types.h> 
#include <sys/stat.h> 
#include <gsl/gsl_rng.h> 
#include <gsl/gsl_randist.h> 
#include <gsl/gsl_permutation.h> 

#include "globals.h" 
#include "mean_field.h" 
#include "net_utils.h" 
#include "con_utils.h" 
#include "rmvnorm.h"

int main(int argc , char** argv) { 
  
  get_args(argc , argv) ;
  
  IF_STRUCTURE = IF_RING || IF_SPEC || IF_LOW_RANK || IF_GAUSS ; 
  
  get_param() ; 
  init_globals() ; 

  create_con_dir() ; 
  
  cout << "coucou" << endl ;

  float *array ; 
  random_normal_multivariate(COV_MAT, array, 4, n_neurons) ;
  
  for(i=0; i<10; i++) { 
    for(int j=0; j<4; j++)
      printf("%f ", array[j+i*4] ) ; 
    printf("\n") ; 
  }
  
  for(i=0; i<n_neurons; i++) { 
    sample[i] = array[0+i*4] ; 
    distractor[i] = array[1+i*4] ; 
    ksi[i] = array[2+i*4] ;
    ksi_1[i] = array[3+i*4] ;      
  }

  cout << "ksi: " ;
  for(i=0; i<10; i++)
    cout << ksi[i] << " " ;
  cout << endl ;

  cout << "ksi_1: " ;
  for(i=0; i<10; i++)
    cout << ksi_1[i] << " " ;
  cout << endl ;

  double cos_ksi ; 
  cos_ksi = cos_array(sample, distractor, n_per_pop[0]) ; 

  cout << " angle_sample_dist: " << acos(cos_ksi) * 180.0 / M_PI << "°" << endl ; 

  cos_ksi = cos_array(ksi, ksi_1, n_per_pop[0]) ; 
  
  // cout << "cos_ksi_ksi_1: " << cos_ksi ; 
  cout << " angle_ksi_ksi_1: " << acos(cos_ksi) * 180.0 / M_PI << "°" << endl ; 
  
  cos_ksi = cos_array(ksi, sample, n_per_pop[0]) ; 
  
  // cout << "cos_ksi_sample: " << cos_ksi ; 
  cout << " angle_ksi_sample: " << acos(cos_ksi) * 180.0 / M_PI << "°" << endl ; 

  cos_ksi = cos_array(ksi, distractor, n_per_pop[0]) ; 

  // cout << "cos_ksi_dist: " << cos_ksi ; 
  cout << " angle_ksi_dist: " << acos(cos_ksi) * 180.0 / M_PI << "°" << endl ; 

  cos_ksi = cos_array(ksi_1, sample, n_per_pop[0]) ; 

  // cout << "cos_ksi_sample: " << cos_ksi ; 
  cout << " angle_ksi1_sample: " << acos(cos_ksi) * 180.0 / M_PI << "°" << endl ; 
  
  cos_ksi = cos_array(ksi_1, distractor, n_per_pop[0]) ; 

  // cout << "cos_ksi_dist: " << cos_ksi ; 
  cout << " angle_ksi1_dist: " << acos(cos_ksi) * 180.0 / M_PI << "°" << endl ; 
  
  // init_ksi() ; 
  
  write_to_file(ksi_path, "ksi", ksi , n_neurons) ; 
  // write_to_file(con_path, "ksi_scaled", ksi_scaled , n_per_pop[0]) ; 
  write_to_file(ksi_path, "sample", sample , n_neurons) ; 
  write_to_file(ksi_path, "dist", distractor , n_neurons) ; 
    
  if(RANK==2) {
    // init_ksi_1() ; 
    write_to_file(ksi_path, "ksi_1", ksi_1 , n_neurons) ; 
    // write_to_file(con_path, "ksi_1_scaled", ksi_1_scaled , n_per_pop[0]) ; 
  }
  
  delete_globals() ; 
  
}
