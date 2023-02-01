#ifndef __PERM__
#define __PERM__

#define FIRST_ITEM 0 ; //firts item in the permu 
long double **count_ ; 
unsigned long n_ ;

void dist_decomp_vector2perm(unsigned long* vec, unsigned long* sigma) {

  unsigned long val, *ident, index ; 
  
  ident = new unsigned long[n_]() ; 
  
  for(i=0; i<n_; i++) 
    ident[i]=i ; 
  
  for(i=0; i<n_-1; i++) {
    val = vec[i] ; 
    index = 0 ;
    
    while( ! ( ident[index]!=-1 && val==0 ) )      
      if(ident[ index ++ ] != -1) 
	val-- ; 
    
    sigma[i] = index + FIRST_ITEM ; 
    ident[index] = -1 ; 
  }
  
  index = 0 ; 
  while(ident[index]==-1)
    index++; 
  
  sigma[n_-1] = index + FIRST_ITEM ; 
  // delete [] ident ; 
}

void random_permu_at_dist_d(unsigned long dist, unsigned long *sigma) { 
  
  double *acum, bound ; 
  unsigned long *v, min, rest_max_dist, pos ; 
  
  acum = new double[n_]() ;  
  v = new unsigned long[n_]() ;
  
  v[n_-1]=0 ;
  
  // generate random distance decomposition vector (inversion vector)
  for(i=0; (i< n_ && dist>0); i++) {
    
    rest_max_dist = (n_ - i - 1 ) * ( n_ - i - 2 ) / 2 ; //con los restantes n' puedes tener distMAx de binom(n`)(2)
	
    if(rest_max_dist >= dist) 
      acum[0] = count_[n_-i-1][dist] ; 
    else
      acum[0] = 0 ; 
    
    min = (n_-i < dist + 1 ) ? (n_ - i ) : dist + 1 ; /////MIN 
    
    for(j=1; j<min; j++ )      
      if(rest_max_dist + j >= dist)
	acum[j] = acum[j-1] + count_[ n_ - i - 1 ] [ dist - j] ; 
      else
	acum[ j ] = 0 ; 
    
    bound = (double)rand() / (double)(RAND_MAX) * acum[min-1] ; 
    pos = 0 ; 
    while(acum[pos] <= bound) 
      pos++;
    
    dist -= pos ;
    v[i] = pos ;
  }
  
  for(j=i; j<n_; j++) 
    v[j] = 0; //the last n-i positions
  
  dist_decomp_vector2perm(v, sigma) ; 
    
  // delete [] v;
  // delete [] acum;
}

void ini_count(){
  // initialize the matrix of the count of permutations at each distance
  //https://oeis.org/A008302

  unsigned long i,j ;
  /* count_ = (long double**) malloc( (unsigned long) (n_+1) * sizeof(long double*) ) ;    */
  count_= new long double*[n_+1]();
  
  for(i=0; i<n_+1; i++)
    count_[i] = new long double[n_ * ( n_- 1 ) / 2 + 1 ]() ;
    /* count_[i] = (long double*) malloc( (unsigned long) (n_ * ( n_- 1 ) / 2 + 1 ) * sizeof(long double) ) ;  */
  
  for(i=0; i<n_+1; i++) {
    count_[i][0] = 1 ;
    
    // for(j=1; j<n_*(n_-1)/2+1; j++) 
    //   count_[i][j] = 0 ; 
  }
  
  for(i=1; i<n_+1; i++){
    for(j=1; j<i*(i-1)/2+1; j++) {
      
      if( j-i >= 0)
	count_[i][j] = count_[i][j-1] + count_[i-1][j] - count_[i-1][j-i];
      else
	count_[i][j] = count_[i][j-1] + count_[i-1][j] ; 
    }
  }
}

void get_perm_at_d(unsigned long N, double dist, unsigned long *sigma) { 

  n_ = N ; 
  unsigned long max_dist = (unsigned long) ( n_*(n_-1)/2 ) ; // maximum distance for two permutations of length n_ 
  unsigned long distance = (unsigned long) ( dist * (double) max_dist ) ; 
  
  /* sigma = new unsigned long[n_]() ;   */
  
  if(distance>max_dist) { 
    cout << "ERROR DIST > MAX_DIST" << endl ; 
    return ; 
  }
  
  cout <<"Random permutation at distance "<< (double) dist << " (" << (double) distance << ")"<< " from the identity " << endl ; 

  ini_count() ; 
  random_permu_at_dist_d(distance, sigma) ; 
  for(j=0; j<100 ; j++)
    cout << sigma[j] << " " ; 
  cout << endl ; 
  
}

#endif 
