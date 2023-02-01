#ifndef __LOWRANKUTILS__
#define __LOWRANKUTILS__

/* #include <gsl/gsl_rng.h>  */
/* #include <gsl/gsl_randist.h>  */
/* #include <gsl/gsl_permutation.h> */

/* #include "rmvnorm.h" */

//////////////////////////////////////////
// Low-rank globals
//////////////////////////////////////////
#define IF_LOW_RANK 0 

const int low_synapse[2] = {1, 0} ; 
const double MEAN_XI[2] = {4.0, 4.0} ; 
const double VAR_XI[2] = {0.5,1.0} ; 

#define IF_LEFT_RIGHT 0 

#define MEAN_XI_LEFT -0.0 
#define VAR_XI_LEFT 5.0 

#define MEAN_XI_RIGHT -0.0 
#define VAR_XI_RIGHT 5.0 

#define RHO 1.0 

#define FIXED_XI 1 
#define SEED_XI 1 

#define IF_FF 0 
#define MEAN_FF 0.0 
#define VAR_FF 0.1 

#define IF_RHO_FF 0  
#define RHO_FF_XI 0.80  
#define VAR_ORTHO 0.0 

#define VAR_FF_LEFT 5.0 
#define VAR_FF_RIGHT 5.0 

#define FIXED_FF 1 

//////////////////////////////////////////
// Low-rank functions 
//////////////////////////////////////////

void create_path_low_rank(string &path, int n_pop) {

  char char_mean_xi[10] ;
  string str_mean_xi ;

  char char_var_xi[10] ;
  string str_var_xi ;
  
  for(int i=0;i<n_pop;i++) {
    sprintf(char_mean_xi, "%0.2f", MEAN_XI[i]) ;
    str_mean_xi = string(char_mean_xi) ;
  
    sprintf(char_var_xi, "%0.2f", VAR_XI[i]) ;
    str_var_xi = string(char_var_xi) ;

    if(n_pop==1)
      path += "/low_rank/xi_mean_" + str_mean_xi + + "_var_" + str_var_xi ;
    else
      if(i==0)
	path += "/low_rank/xi_E_mean_" + str_mean_xi + + "_var_" + str_var_xi ;
      else
	path += "_xi_I_mean_" + str_mean_xi + + "_var_" + str_var_xi ;
  }
} 

void create_path_low_rank_left_right(string &path) {

  char char_mean_xi_left[10] ;
  sprintf(char_mean_xi_left, "%0.2f", MEAN_XI_LEFT) ;
  string str_mean_xi_left = string(char_mean_xi_left) ;
  
  char char_var_xi_left[10] ;
  sprintf(char_var_xi_left, "%0.2f", VAR_XI_LEFT) ;
  string str_var_xi_left = string(char_var_xi_left) ;
  
  char char_mean_xi_right[10] ;
  sprintf(char_mean_xi_right, "%0.2f", MEAN_XI_RIGHT) ;
  string str_mean_xi_right = string(char_mean_xi_right) ;
  
  char char_var_xi_right[10] ;
  sprintf(char_var_xi_right, "%0.2f", VAR_XI_RIGHT) ;
  string str_var_xi_right = string(char_var_xi_right) ;

  char char_rho[10] ;
  sprintf(char_rho, "%0.2f", RHO) ;
  string str_rho = string(char_rho) ;

  path += "/low_rank/xi_left_mean_" + str_mean_xi_left + + "_var_" + str_var_xi_left ;
  path += "_xi_right_mean_" + str_mean_xi_right + + "_var_" + str_var_xi_right ;
  path += "/rho_"+ str_rho ;
}

void create_path_ff(string &path) {

  char char_mean_ff[10] ;
  sprintf(char_mean_ff, "%0.2f", MEAN_FF) ;
  string str_mean_ff = string(char_mean_ff) ;
  
  char char_var_ff[10] ;
  sprintf(char_var_ff, "%0.2f", VAR_FF) ;
  string str_var_ff = string(char_var_ff) ;
  
  char char_var_ortho[10] ;
  
  if(IF_RHO_FF) 
    sprintf(char_var_ortho, "%0.2f", RHO_FF_XI) ; 
  else
    sprintf(char_var_ortho, "%0.2f", VAR_ORTHO) ; 	 
  
  string str_var_ortho = string(char_var_ortho) ; 

  if(IF_RHO_FF)
    path += "/ff_mean_" + str_mean_ff + "_var_" + str_var_ff + "_rho_" + str_var_ortho ; 
  else
    path += "/ff_mean_" + str_mean_ff + "_var_" + str_var_ff + "_ortho_" + str_var_ortho ; 
}

void create_path_ff_left_right(string &path) {
  char char_mean_ff[10] ;
  sprintf(char_mean_ff, "%0.2f", MEAN_FF) ;
  string str_mean_ff = string(char_mean_ff) ;
  
  char char_var_ff_left[10] ;
  sprintf(char_var_ff_left, "%0.2f", VAR_FF_LEFT) ;
  string str_var_ff_left = string(char_var_ff_left) ;
  
  char char_var_ff_right[10] ;
  sprintf(char_var_ff_right, "%0.2f", VAR_FF_RIGHT) ;
  string str_var_ff_right = string(char_var_ff_right) ;
  
  char char_var_ortho[10] ;
  sprintf(char_var_ortho, "%0.2f", VAR_ORTHO) ;	 
  string str_var_ortho = string(char_var_ortho) ;
  
  path += "/ff_mean_" + str_mean_ff + "_var_left_" + str_var_ff_left +"_right_" + str_var_ff_right  + "_ortho_" + str_var_ortho ; 
}

void multivariate_ksi_stimuli(int size) {
  time_t initime ;

  // alloc gsl vectors and cov matrix 
  gsl_vector *x = gsl_vector_calloc(size),
    *mean = gsl_vector_calloc(size);
  gsl_matrix *m = gsl_matrix_alloc(size, size),
    *rm = gsl_matrix_alloc(size, size);
  gsl_rng *r;

  // init gsl seed 
  r = gsl_rng_alloc(gsl_rng_mt19937) ; 
  gsl_rng_set(r, initime) ; 
  printf("The SEED is %li.\n\n",initime) ; 
  
  // init gsl vectors and cov matrix 
  for(int i=0; i<size; i++) {
    gsl_vector_set(x, i, 0.0) ; 
    for(int j=0; j<size; j++) 
      gsl_matrix_set(m, i, j, COV_MAT[i+j*size]) ; 
  }

  FILE *f;
  f = fopen(path + "/gsl_rmvn_vec.txt", "w") ; 
  
  for(int k=0; k<n_neurons; k++){ 
    rmvnorm(r, size, mean, m, x) ; 
    if(k<10) 
      for(int i=0; i<size; i++) 
	cout << gsl_vector_get(x, i) << " " ; 
    cout << endl ;
    
    for(int i=0; i<size; i++) 
      fprintf(f, "%g ", gsl_vector_get(x,i) ) ; 
    fprintf(f, "\n") ; 
    
  } 
  
} 

#endif
