#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
/* ----------------------------- */
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
/* ----------------------------- */

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

float mean_array(float *array, int size) {
  float mean = 0 ; 
  for(int i=0; i<size; i++)
    mean += array[i] ;
  mean /= size ;

  return mean ;
}

float var_array(float *array, int size) {
  float mean = 0 ; 
  float var = 0 ;
  for(int i=0; i<size; i++)
    var += array[i]*array[i] ; 
  var /= size ; 
  
  mean = mean_array(array, size) ; 
  var -= mean * mean ;

  return var ; 
}

float covar_array(float *a, float *b, int size) {
  float ma=0, mb = 0 ; 
  float covar = 0 ;
  
  for(int i=0; i<size; i++)
    covar += a[i]*b[i] ; 
  covar /= size ; 
  
  ma = mean_array(a, size) ;  
  mb = mean_array(b, size) ;
  
  covar -= ma * mb ;
  
  return covar ; 
}

float dmvnorm(const int n, const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var){
  /* multivariate normal density function    */
  /*
   *    n       dimension of the random vetor
   *    mean    vector of means of size n
   *    var     variance matrix of dimension n x n
   */
  int s;
  double ax,ay; // must be double 
  gsl_vector *ym, *xm;
  gsl_matrix *work = gsl_matrix_alloc(n,n),
    *winv = gsl_matrix_alloc(n,n);
  gsl_permutation *p = gsl_permutation_alloc(n);
 
  gsl_matrix_memcpy( work, var );
  gsl_linalg_LU_decomp( work, p, &s );
  gsl_linalg_LU_invert( work, p, winv );
  ax = gsl_linalg_LU_det( work, s );
  gsl_matrix_free( work );
  gsl_permutation_free( p );
 
  xm = gsl_vector_alloc(n);
  gsl_vector_memcpy( xm, x);
  gsl_vector_sub( xm, mean );
  ym = gsl_vector_alloc(n);
  gsl_blas_dsymv(CblasUpper,1.0,winv,xm,0.0,ym);
  gsl_matrix_free( winv );
  gsl_blas_ddot( xm, ym, &ay);
  gsl_vector_free(xm);
  gsl_vector_free(ym);
  ay = exp(-0.5*ay)/sqrt( pow((2.0*M_PI),n)*ax );
 
  return ay;
}


int rmvnorm(const gsl_rng *r, const int n, const gsl_vector *mean, const gsl_matrix *var, gsl_vector *result){
  /* multivariate normal distribution random number generator */
  /*
   *    n       dimension of the random vetor
   *    mean    vector of means of size n
   *    var     variance matrix of dimension n x n
   *    result  output variable with a sigle random vector normal distribution generation
   */
  int k;
  gsl_matrix *work = gsl_matrix_alloc(n,n);
 
  gsl_matrix_memcpy(work,var);
  gsl_linalg_cholesky_decomp(work);
 
  for(k=0; k<n; k++)
    gsl_vector_set( result, k, gsl_ran_ugaussian(r) );
  
  gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, result);
  gsl_vector_add(result,mean);
 
  gsl_matrix_free(work);
 
  return 0;
}


void random_normal_multivariate(float *cov_mat, float* &array, int size, int items) {
  
  time_t initime ;
  
  // alloc gsl vectors and cov matrix 
  gsl_vector *x = gsl_vector_calloc(size),
    *mean = gsl_vector_calloc(size);
  gsl_matrix *m = gsl_matrix_alloc(size, size) ;
  gsl_rng *r;
  
  // init gsl seed 
  r = gsl_rng_alloc(gsl_rng_mt19937) ;
  if(FIX_KSI_SEED){
    gsl_rng_set(r, (long int) SEED_KSI) ; 
    printf("The SEED is %li.\n", (long int) SEED_KSI) ;
  }
  else {
    gsl_rng_set(r, initime) ; 
    printf("The SEED is %li.\n", initime) ; 
  }
  // init gsl vectors and cov matrix 
  for(int i=0; i<size; i++) {
    gsl_vector_set(x, i, 0.0) ; 
    gsl_vector_set(mean, i, 0.0) ; 
    for(int j=0; j<size; j++) 
      gsl_matrix_set(m, i, j, cov_mat[i+j*size]) ; 
  }
  
  FILE *f;
  f = fopen("gsl_rmvn_vec.txt", "w") ; 
  
  for(int k=0; k<items; k++){ 
    rmvnorm(r, size, mean, m, x) ; 
    
    for(int i=0; i<size; i++) { 
      /* fprintf(stdout, "%f ", gsl_vector_get(x,i) ) ;  */
      fprintf(f, "%f ", gsl_vector_get(x,i) ) ; 
    }
    /* fprintf(stdout, "\n") ;  */
    fprintf(f, "\n") ; 
  }
  
  fclose(f);
  
  // float array[size*items] ;
  int dum ;
  array = (float *) malloc( (int) (size * items) * sizeof(float) ) ;
  
  f = fopen("gsl_rmvn_vec.txt", "r") ;
  // dum = fread(&array, sizeof array, 1, file) ;
  int n=0 ;
  while(fscanf(f, "%f", &array[n])!=EOF) n++ ;
  fclose(f) ;
  
} 
