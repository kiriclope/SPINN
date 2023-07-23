#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <sys/stat.h>
#include <libconfig.h++>

#include "globals.h"
#include "mean_field.h"
#include "net_utils.h"

#include "con_utils.h"
#include "con_LR_utils.h"
#include "con_ff_utils.h"

#include "mat_utils.h"

#include "stp_utils.h"
#include "tasks_utils.h"

#include "lif_utils.h"
// #include "binary_utils.h"


int main(int argc , char** argv) {

  get_args(argc , argv) ;

  IF_STRUCTURE = IF_RING || IF_SPEC || IF_LOW_RANK || IF_GAUSS ;
  int IF_SIM= IF_BIN || IF_LIF || IF_RATE ;

  IF_STIM = IF_DPA || IF_DUAL || IF_DRT || IF_STEP || IF_CHRISTOS ;
  IF_NOISE = IF_COS_NOISE || IF_GAUSS_NOISE ;

  cout << "STIM" << IF_STIM << " IF_NOISE " << IF_NOISE << endl ;

  read_params() ;
  get_param() ;

  init_globals() ;

  if(IF_LOW_RANK)
    if(IF_GEN_KSI)
      gen_ksi() ;

  if(IF_GEN_CON)
    gen_con_sparse_vec() ;
  else
    if(IF_LIF)
      get_con_sparse_vec() ;

  if(IF_FF_LAYER) {
    if(IF_SPARSE_FF)
      if(IF_GEN_CON_FF)
        gen_con_sparse_vec_ff() ;
      else
        get_con_sparse_vec_ff() ;
    else
      init_theta_ff() ;
  }

  if(IF_SIM) {

    create_dir() ;
    mean_field_rates() ;

    if(IF_LIF)
      run_sim_lif() ;

    // if(IF_BIN)
    //   run_sim_bin() ;
  }
}
