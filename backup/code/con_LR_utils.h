#ifndef _LR_UTILS_ 
#define _LR_UTILS_

#include <gsl/gsl_rng.h> 
#include <gsl/gsl_randist.h> 
#include <gsl/gsl_permutation.h>

#include "rmvnorm.h"

/* #include "low_rank_utils.h" */

float norm_array(float *array, size_t n) { 
  
  float norm = 0 ; 
  for(i=0; i<n; i++) 
    norm += array[i]*array[i] ;
  
  return sqrt(norm) ;
}

float cos_array(float *a, float *b, size_t n) {
  float cosine_ab = 0, norm_a=0, norm_b=0 ; 
  for(i=0; i<n; i++) 
    cosine_ab += a[i]*b[i] ; 

  norm_a = norm_array(a, n) ;
  norm_b = norm_array(b, n) ;

  return cosine_ab / norm_a / norm_b ; 
}

void outer_product(float *a, float *b, size_t n, float *&outer) { 
  
  outer = new float [n*n]() ; 
  for(i=0; i<n; i++)
    for(j=0; j<n; j++)
      outer[i+j*n] = a[i]*b[j] ; 
  
}

void normalize_array(float *&a, size_t n) {
  float norm = 0 ;
  norm = norm_array(a,n) ;
  
  for(i=0; i<n; i++)
    a[i] /= norm ; 
}
  
void rotate_ab(float* a, float* b, size_t n, float angle, float *&result) {

  normalize_array(a, n) ; 
  normalize_array(b, n) ; 
  
  float *Id;
  Id = new float [n*n]() ; 
  for(i=0; i<n; i++)
    for(j=0; j<n; j++)
      Id[i+j*n] = 1.0 ; // Identity Matrix 

  float *outer_ab, *outer_ba, *outer_aa, *outer_bb; 
  outer_product(a, a, n, outer_aa) ;
  outer_product(b, b, n, outer_bb) ;  
  outer_product(a, b, n, outer_ab) ;
  outer_product(b, a, n, outer_ba) ;
  
  /* R = I + ( np.outer(n2,n1) - np.outer(n1,n2) ) * np.sin(a) + ( np.outer(n1,n1) + np.outer(n2,n2) ) * (np.cos(a)-1) ; */
  
  float *R;
  R = new float [n*n]() ;
  
  for(i=0; i<n; i++)
    for(j=0; j<n; j++)
      R[i+j*n] = Id[i+j*n] + ( outer_ba[i+j*n] - outer_ab[i+j*n] ) * sin(angle) + ( outer_aa[i+j*n] + outer_bb[i+j*n] ) * (cos(angle) - 1.0) ;

  for(i=0; i<n; i++)
    for(j=0; j<n; j++)
      result[i] += R[i+j*n] * b[j] ; 
  
  delete [] Id ;
  delete [] R ; 
}

void angle_ksi() { 
   
  float cos_ksi_0 = 0.0 ;
  cos_ksi_0 = cos_array(ksi, ksi_1, n_per_pop[0]) ; 
  
  cout << "cos_ksi_0: " << cos_ksi_0 ; 
  cout << " angle_ksi_0: " << acos(cos_ksi_0) * 180.0 / M_PI << "°" << endl ; 
  
  /* float norm_ksi=0, norm_ksi_1=0, dot=0 ;  */
  /* norm_ksi = norm_array(ksi, n_per_pop[0]) ;  */
  /* norm_ksi_1 = norm_array(ksi, n_per_pop[0]) ;  */
  
  /* for(i=0; i<n_per_pop[0]; i++)  */
  /*   dot += ksi[i]*ksi_1[i] ;  */
  
  /* for(i=0; i<n_per_pop[0]; i++)  */
  /*   ksi_1[i] = ksi_1[i] - dot/norm_ksi/norm_ksi * ksi[i] ; // Gram Schmidt ortho  */
  /*   /\* ksi_1[i] = ksi_1[i] - dot/norm_ksi_1/norm_ksi_1 * ksi_1[i] ; // Gram Schmidt ortho  *\/  */
  
  /* float cos_ksi = 0.0 ;  */
  /* cos_ksi = cos_array(ksi, ksi_1, n_per_pop[0] ) ;  */
  
  /* cout << "cos_ksi: " << cos_ksi ;  */
  /* cout << " angle_ksi: " << acos(cos_ksi) * 180.0 / M_PI << "°" << endl ; */
  
  /* float *ksi_rotate; */
  /* ksi_rotate = new float [n_per_pop[0]]() ;  */
  
  /* rotate_ab( ksi, ksi_1, n_per_pop[0], ANGLE_KSI, ksi_rotate) ;  */
  
  /* cos_ksi = cos_array(ksi, ksi_rotate, n_per_pop[0] ) ;  */
  
  /* cout << "cos_ksi: " << cos_ksi ;  */
  /* cout << " angle_ksi: " << acos(cos_ksi) * 180.0 / M_PI << "°" << endl ; */
  
}

void angle_maps() { 

  cout << "theta: " ; 
  for(i=0; i<10; i++) 
    cout << theta[i] << " " ; 
  cout << endl ;
  
  cout << "theta_1: " ; 
  for(i=0; i<10; i++) 
    cout << theta_1[i] << " " ; 
  cout << endl ; 
  
  float sum_theta=0.0, norm_theta ; 
  
  for(i=0; i<n_per_pop[0]; i++) 
    sum_theta += theta[i] ; 
  
  cout << "sum theta: " << sum_theta << endl ;
  float sum_theta_1=0.0, norm_theta_1 ;
  
  for(i=0; i<n_per_pop[0]; i++) { 
    theta_1[i] = theta_1[i] - theta[i] * theta_1[i] / sum_theta ; 
    theta_1[i] = theta_1[i] + cos(MAP_ANGLE) / sum_theta ; 
    /* sum_theta_1 += theta_1[i] ;  */ 
  }

  cout << "theta_1 ortho: " ; 
  for(i=0; i<10; i++) 
    cout << theta_1[i] << " " ; 
  cout << endl ; 
  
  float cos_maps = 0.0 ; 
  for(i=0; i<n_per_pop[0]; i++) {
    /* theta_1[i] += cos(MAP_ANGLE) / sum_theta_1 ;  */ 
    cos_maps += theta[i]*theta_1[i] ; 
  }
  
  norm_theta = norm_array(theta, n_per_pop[0]) ;
  norm_theta_1 = norm_array(theta_1, n_per_pop[0]) ; 
  
  cos_maps = cos_maps / norm_theta / norm_theta_1 ; 
  
  cout<< "norm_theta " << norm_theta << " norm_theta_1 " << norm_theta_1 << endl ; 
  cout << "cos_maps: " << cos_maps << ", angle_maps: " << acos(cos_maps) * 180.0 / M_PI << "°" << endl ;
}

void my_shuffle(void *base, size_t n, size_t size) { 
  
  const gsl_rng_type * T ; 
  gsl_rng * r ; 
  T = gsl_rng_default ; 
  r = gsl_rng_alloc (T) ;
  
  if(!FIX_MAP_SEED)
    gsl_rng_set(r, clock());
  else {
    cout << "MAP SEED " << exp( (float) MAP_SEED) << endl ; 
    gsl_rng_set(r, exp( (float) MAP_SEED) ); 
  }
  gsl_ran_shuffle (r, base, (size_t) n, size) ; 
  
}

size_t* permutation(size_t N) {

  FILE *pFile ; 
  string file_path ; 
  
  /* file_path = con_path + "idx_perm.dat" ;  */
  /* cout << file_path << endl ;  */
  
  /* pFile = fopen(file_path.c_str(), "wb") ;  */
  
  const gsl_rng_type * T ; 
  gsl_rng * r ; 
  gsl_permutation * p = gsl_permutation_alloc (N) ; 
  
  gsl_rng_env_setup(); 
  T = gsl_rng_default; 
  r = gsl_rng_alloc (T); 
  gsl_rng_set(r, clock()); 
  
  gsl_permutation_init (p); 
  gsl_ran_shuffle (r, p->data, N, sizeof(size_t)) ; 
  /* gsl_permutation_fprintf (pFile, p, " %u");  */
  
  /* fclose(pFile) ;  */
  return p->data ; 
}

float distance_between_maps() {
  float sum = 0 ;
  for(i=0; i<n_per_pop[0]; i++) 
    if( fmod(2.0 * M_PI * (float) (idx_perm_E[i]-i) / (float) n_per_pop[0], 2.0*M_PI) >= M_PI/4.0 ) 
      sum += 1 ; 
  
  return sum / (float) n_per_pop[0] ; 
}


void init_theta_1() {
  
  /* idx_perm = permutation( (size_t) n_per_pop[0] ) ; */
  
  for(i=0; i<n_per_pop[0]; i++) 
    idx_perm_E[i] = i ; 
  
  my_shuffle(idx_perm_E, n_per_pop[0], sizeof(unsigned long) ) ; 
  
  /* float distance = 0 ; */
  /* distance = distance_between_maps() ; */
  /* cout << "distance between maps " << distance << endl ; */
  
  cout << "idx_perm_E: " ; 
  for(i=0;i<10;i++) 
    cout << idx_perm_E[i] << " " ; 
  cout << endl ; 
  
  /* for(i=0; i<n_per_pop[1]; i++)  */
  /*   idx_perm_I[i] = i + n_per_pop[0] ;  */
  
  /* my_shuffle(idx_perm_I, n_per_pop[1], sizeof(unsigned long) ) ;  */
  
  /* cout << "idx_perm_I: " ;  */
  /* for(i=0;i<10;i++)  */
  /*   cout << idx_perm_I[i] << " " ;  */
  /* cout << endl ;  */
  
  for(i=0; i<n_pop; i++)
    if(i==0)
      for(j=cum_n_per_pop[i]; j<cum_n_per_pop[i+1]; j++) 
	idx_perm[j] = idx_perm_E[j] ; 
  /* else  */
  /*   for(j=cum_n_per_pop[i]; j<cum_n_per_pop[i+1]; j++)  */
  /* 	idx_perm[j] = idx_perm_I[j-cum_n_per_pop[i]] ;  */
  
  cout << "idx_perm: " ; 
  for(i=0;i<10;i++) 
    cout << idx_perm[i] << " " ; 
  cout << endl ; 
  
  for(i=0; i<n_pop; i++) 
    for(j=cum_n_per_pop[i]; j<cum_n_per_pop[i+1]; j++)
      theta_1[j] = theta[idx_perm[j]] ; 
  
  /* if(i==0)  */
  /* 	theta_1[j] = 2.0 * M_PI * (float) (idx_perm[j]-cum_n_per_pop[i]) / (float) n_per_pop[i] ;  */
  /* else  */
  /* 	theta_1[j] = 2.0 * M_PI * (float) (j-cum_n_per_pop[i]) / (float) n_per_pop[i] ;  */
  
  /* angle_maps() ;  */


  cout << "theta_1: " ; 
  for(i=0; i<10; i++) 
    cout << theta_1[i] << " " ; 
  cout << endl ;   
}

void create_LR_con_dir() {
  
  ostringstream str_seed_ksi ; 
  str_seed_ksi << fixed << setprecision(0) << SEED_KSI ; 
  
  ksi_path += "../connectivity/" + to_string(n_pop) +"pop" ;
  ksi_path += "/NE_" + to_string(n_per_pop[0]/1000) +  "_NI_" + to_string(n_per_pop[1]/1000) ; 
  ksi_path += "/low_rank/rank_" + to_string(RANK) ; 
  
  if(FIX_KSI_SEED) 
    ksi_path += "/seed_ksi_" + str_seed_ksi.str() ; 
  
  make_dir(ksi_path) ; 
} 

void gen_ksi() {

  cout << "Generate ksi " << endl ; 
  create_LR_con_dir() ; 
  
  float *array ; 
  random_normal_multivariate(COV_MAT, array, 4, n_per_pop[0]) ; 
  
  /* for(i=0; i<10; i++) {  */
  /*   for(int j=0; j<4; j++)  */
  /*     printf("%f ", array[j+i*4] ) ;  */
  /*   printf("\n") ;  */
  /* } */
  
  for(i=0; i<n_per_pop[0]; i++) { 
    sample_A[i] = array[0+i*4] ;
    sample_B[i] = array[1+i*4] ;
    ksi[i] = array[2+i*4] ;
    if(RANK==2)
      ksi_1[i] = array[3+i*4] ;
  }
  
  cout << "sample A: " ;
  cout << "mean " << mean_array(sample_A, n_per_pop[0]) << " var " << var_array(sample_A, n_per_pop[0]) << endl ; 
  
  cout << "sample B: " ;
  cout << "mean " << mean_array(sample_B, n_per_pop[0]) << " var " << var_array(sample_B, n_per_pop[0]) << endl ; 

  cout << "covar sample A/B " << covar_array(sample_A, sample_B, n_per_pop[0]) << endl ; 
  
  cout << "ksi: " ;
  cout << "mean " << mean_array(ksi, n_per_pop[0]) << " var " << var_array(ksi, n_per_pop[0]) << endl ; 
  cout << "covar ksi/samples " << covar_array(ksi, sample_A, n_per_pop[0]) << " | " << covar_array(ksi, sample_B, n_per_pop[0]) << endl ; 
  
  if(RANK==2) {
    cout << "ksi_1: " ;
    cout << "mean " << mean_array(ksi_1, n_per_pop[0]) << " var " << var_array(ksi_1, n_per_pop[0]) << endl ;
    cout << "covar ksi/samples " << covar_array(ksi_1, sample_A, n_per_pop[0]) << " | " << covar_array(ksi_1, sample_B, n_per_pop[0]) << endl ; 
    cout << "covar ksi/ksi_1 " << covar_array(ksi, ksi_1, n_per_pop[0]) << endl ;
    
    angle_ksi() ; 
  }
    
  cout << "###############################################" << endl ;
  
  write_to_file(ksi_path, "ksi", ksi , n_per_pop[0]) ;
  
  if(RANK==2)
    write_to_file(ksi_path, "ksi_1", ksi_1 , n_per_pop[0]) ; 
  
  write_to_file(ksi_path, "sample_A", sample_A , n_per_pop[0]) ; 
  write_to_file(ksi_path, "sample_B", sample_B , n_per_pop[0]) ; 
  
}


void get_ksi(){
  create_LR_con_dir() ;
  
  read_from_file(ksi_path, "ksi", ksi, n_per_pop[0]) ; 
  
  read_from_file(ksi_path, "sample_A", sample_A, n_per_pop[0]) ; 
  read_from_file(ksi_path, "sample_B", sample_B, n_per_pop[0]) ; 
  
  cout << "ksi " << endl ;
  for(i=0;i<10;i++)
    cout << ksi[i] << " " ;
  cout << endl ;
  
  cout << "sample_A " << endl ;
  for(i=0;i<10;i++)
    cout << sample_A[i] << " " ;
  cout << endl ;

  cout << "sample_B " << endl ;
  for(i=0;i<10;i++)
    cout << sample_B[i] << " " ;
  cout << endl ;

  if(RANK==2) {
    read_from_file(ksi_path, "ksi_1", ksi_1, n_per_pop[0]) ;
    
    cout << "ksi_1 " << endl ;
    for(i=0;i<10;i++)
      cout << ksi_1[i] << " " ;
    cout << endl ;
  }  
  
  if(IF_SPEC && RANK==2) {
    read_from_file(con_path, "idx_perm", idx_perm, n_neurons) ; 
    read_from_file(con_path, "theta_1", theta_1, n_neurons) ; 
    
    cout << "idx_perm: " ; 
    for(i=0;i<5;i++) 
      cout << idx_perm[i] << " " ; 
    cout << endl ;
    
    cout << "theta_1: " ; 
    for(i=0;i<5;i++)
      cout << theta_1[i] << " " ; 
    cout << endl ;    
  } 
  
}

#endif
