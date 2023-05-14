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
  
  gen_con_sparse_vec() ; 
}
