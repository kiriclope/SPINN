#ifndef _FF_UTILS_ 
#define _FF_UTILS_

///////////////////////////////////////////////////////////////////////////
// feed_forward connectivity
///////////////////////////////////////////////////////////////////////////

void init_theta_ff() {
  
  theta_ff = new float [N_POISSON]() ; 
  
  for(i=0; i<N_POISSON; i++) {
    theta_ff[i] = ( 2.0 * M_PI * (float) (i) ) / (float) N_POISSON ;
    /* theta_ff[i] = 2.0 * M_PI *  unif(con_gen) ; */
  }
  
  cout << "theta_ff: " ; 
  for(i=0; i<10; i++) 
    cout << theta_ff[i] << " " ; 
  cout << endl ; 
  
}

void make_con_prob_ff() { 

  con_prob_ff = new float [n_neurons*N_POISSON]() ; 
  K_over_N_FF = K_FF / (float) N_POISSON ; 
  
  if(IF_RING_FF || IF_SPEC_FF || IF_GAUSS_FF || IF_GAUSS_SPEC_FF) { 

    cout << "init theta" << endl ;   
    init_theta() ;
    cout << "init theta_ff" << endl ; 
    init_theta_ff() ;
    
    if(IF_GAUSS_FF) { 
      cout << "gaussian structure: " ; 
      cout <<"sigma FF "<< KAPPA_FF << endl ; 
      
      prefactor_ff = new float [n_neurons*n_pop]() ; 
      
      for(i=0; i<n_neurons; i++) { // need loop over post before pre to compute prefactor
	post_pop = which_pop[i] ;
	
	for(j=0; j<N_POISSON; j++) { // Jij: j (pre) to i (post) 
	    con_prob_ff[j + i * N_POISSON] = Gaussian1D(theta[i]-theta_ff[j], KAPPA_FF) ; 
	    // Pij = Zb Cij with Zb = K / sum_j Cij then sum_j Pij = K 
	    prefactor_ff[post_pop + i * n_pop] += con_prob_ff[j + i * N_POISSON] ; // sum over presynaptic neurons j
	  }
	
	for(j=0; j<n_neurons; j++) 
	  con_prob_ff[j + i * N_POISSON] *= K_FF / prefactor_ff[post_pop + i * n_pop] ; 
      } // end loop post
      
      delete [] prefactor_ff ; 
      
    }// endif Gauss
  
    if(IF_GAUSS_SPEC_FF) { 
      cout << "gaussian spec structure: " ; 
      cout <<"sigma FF "<< KAPPA_FF << endl ; 
      
      prefactor_ff = new float [n_neurons*n_pop]() ; 
      K_over_N_FF = (K_FF-sqrt(K_FF)) / (float) N_POISSON ; 
      
      for(i=0; i<n_neurons; i++) { // need loop over post before pre to compute prefactor
	post_pop = which_pop[i] ;
	
	/* if(post_pop==0) { */
	  for(j=0; j<N_POISSON; j++) { // Jij: j (pre) to i (post) 
	    con_prob_ff[j + i * N_POISSON] = Gaussian1D(theta[i]-theta_ff[j], KAPPA_FF) ; 
	    // Pij = Zb Cij with Zb = K / sum_j Cij then sum_j Pij = K 
	    prefactor_ff[post_pop + i * n_pop] += con_prob_ff[j + i * N_POISSON] ; // sum over presynaptic neurons j 
	  }
	  
	  for(j=0; j<n_neurons; j++) {
	    con_prob_ff[j + i * N_POISSON] *= sqrt(K_FF) / prefactor_ff[post_pop + i * n_pop] ; 
	    con_prob_ff[j + i * N_POISSON] += K_over_N_FF ; 
	  }
	  
	  /* } */
	  /* else */
	  /*   con_prob_ff[j + i * N_POISSON] = K_over_N_FF ;  */
      } // end loop post
      
      delete [] prefactor_ff ; 
      
    }// endif Gauss
  
    if(IF_RING_FF || IF_SPEC_FF) {
      
      if(IF_RING_FF) {
	cout << "ring: " ; 
	cout << "kappa_FF, " << KAPPA_FF << endl ; 
      }
      
      if(IF_SPEC_FF) {
	cout << "spec: " ; 
	cout << "kappa_FF, " << KAPPA_FF << endl ;
      }
    
      float kappa_ff ;
      
      if(IF_RING_FF) 
	kappa_ff = 2.0 * KAPPA_FF * K_over_N_FF ; 
      if(IF_SPEC_FF) 
	kappa_ff = K_over_N_FF * KAPPA_FF / sqrt(K_FF) ; 
      
      if(IF_RING_FF && KAPPA_FF > 0.5) 
	cout << "ERROR: KAPPA_FF TOO LARGE" << endl ; 
      if(IF_SPEC_FF && KAPPA_FF > sqrt(K_FF)) 
	cout << "ERROR: KAPPA_FF TOO LARGE" << endl ; 
      
      /* double Delta ; */
      
      for(j=0; j<N_POISSON; j++) { 
	/* Delta = unif(rand_gen) * 2.0 * M_PI ; */
	for(i=0; i<n_neurons; i++) { 
	  post_pop = which_pop[i] ; 
	  /* kappa_ff = sigma_FF[post_pop] ;  */ 
	  
	  if(i==n_per_pop[0]) {
	    if(IF_SPEC_FF) 
	      kappa_ff = KAPPA_FF / sqrt(K_FF) * K_over_N_FF  ; 
	    if(IF_RING_FF) 
	      kappa_ff = 2.0 * KAPPA_FF * K_over_N_FF ; 
	  }
	  
	  con_prob_ff[j + i * N_POISSON] += K_over_N_FF ; 
	  /* if(i<n_per_pop[0])  */
	  con_prob_ff[j + i * N_POISSON] += kappa_ff * cos( theta[i] - theta_ff[j]) ; 
	  /* con_prob_ff[j + i * N_POISSON] += K_over_N_FF * kappa_ff * sin( theta[i] - theta_ff[j] ) ;  */
	  
	}
      }
    } // end ring
  }
  else {
    cout << "random FF" << endl ; 
    
    for(j=0; j<N_POISSON; j++) { // Pij -> j (cols, pre) to i (rows, post) 
      pre_pop = which_pop[j] ; 
      for(i=0; i<n_neurons; i++) 
	con_prob_ff[j + i * N_POISSON] = K_over_N_FF ; // Pij -> P[j + i * n_cols] (index = indexX * arrayWidth + indexY) ; 
    } 
  } 
}

void make_con_vec_ff() {
  
  cout << "generate FF connectivity vector" << endl ; 
  
  con_vec_ff = new int [N_POISSON*n_neurons]() ; 
  
  for(i=0; i<n_neurons; i++) { 
    for(j=0; j<N_POISSON; j++) { // Jij -> j (cols, pre) to i (rows, post) 
      if(con_prob_ff[j + i * N_POISSON] >= unif(con_gen) ) { 
  	con_vec_ff[j + i * N_POISSON] = 1 ; 
  	total_n_post_ff++ ; 
      } 
    } 
  } 
  
  delete [] con_prob_ff ; 
  
  if(IF_SAVE_CON_VEC && !IF_TRIALS) { 
    cout << "saving connectivity vector to: "; 
    write_to_file(con_path_ff, "con_vec_ff", con_vec_ff , n_neurons*N_POISSON) ; 
  } 
} 


void make_con_sparse_rep_ff() { 
  
  unsigned long counter = 0 ;  
  
  cout << "generate sparse representation FF" << endl ;
  cout << "total_n_post_ff " << total_n_post_ff << endl ; 
  
  idx_post_ff = new unsigned long [N_POISSON]() ; 
  id_post_ff = (unsigned long *) malloc( (unsigned long) total_n_post_ff * sizeof(unsigned long) ) ; 
  
  n_post_ff = new int [N_POISSON]() ;
  
  avg_n_post_ff = new int [n_pop]() ;
  
  for(j=0; j<N_POISSON; j++) { // Jij -> j (cols, pre) to i (rows, post) 
    
    for(i=0; i<n_neurons; i++) { 
      post_pop = which_pop[i] ;
      
      if(con_vec_ff[j + i * N_POISSON]) { // j->i = 1 with proba Kj/Nj 
	id_post_ff[counter] = i ; 
	n_post_ff[j]++ ; // Kb to post pop a, K/N * N = K 
	avg_n_post_ff[post_pop]++ ; 
	counter++ ; 
      } 
    } 
  } 
  
  delete [] con_vec_ff ; 
  
  cout << "average number of postsynaptic connections per FF: " ;
  for(i=0; i<n_pop; i++)
    cout << i << " " << avg_n_post_ff[i] / n_per_pop[i] << " " ;
  cout << endl ;
  
  delete [] avg_n_post_ff ;
  
  cout << "total number of postsynaptic connections per FF neuron: " << endl ;
  
  for(k=N_POISSON; k>N_POISSON-10; k--)
    cout << n_post_ff[k] << " " ; 
  cout << endl ;
  
  
  idx_post_ff[0] = 0 ; 
  for(i=1; i<N_POISSON; i++) 
    idx_post_ff[i] = idx_post_ff[i-1] + n_post_ff[i-1] ; 
  
}

void create_ff_con_dir() {
  
  con_path_ff += "/NE_" + to_string(n_per_pop[0]/1000) +  "_NI_" + to_string(n_per_pop[1]/1000) ;   
  con_path_ff += "/N_FF_" + to_string(N_POISSON/1000) ;   
  con_path_ff += "/K_FF_" + to_string( (int) K_FF) ;
  
  ostringstream str_kappa_ff ;
  str_kappa_ff << fixed << setprecision(2) << KAPPA_FF ;
  if(IF_RING_FF) 
    con_path_ff += "/ring/kappa_FF_" + str_kappa_ff.str() ; 
  if(IF_SPEC_FF) 
    con_path_ff += "/spec/kappa_FF_" + str_kappa_ff.str() ; 
  if(IF_GAUSS_FF) 
    con_path_ff += "/gauss/sigma_FF_" + str_kappa_ff.str() ;
  if(IF_GAUSS_SPEC_FF) 
    con_path_ff += "/gauss_spec/sigma_FF_" + str_kappa_ff.str() ;    
  
  make_dir(con_path_ff) ; 
  cout << "created directory: " ; 
  cout << con_path_ff << endl ; 
  
}

void save_ff_to_con_file() {

  cout << "save FF to file" << endl ;
    
  write_to_file(con_path_ff, "id_post_ff", id_post_ff, total_n_post_ff) ; 
  write_to_file(con_path_ff, "idx_post_ff", idx_post_ff, N_POISSON) ; 
  write_to_file(con_path_ff, "n_post_ff", n_post_ff, N_POISSON) ; 
  
}

void gen_con_sparse_vec_ff() {
  
  cout << "ff connectivity path" << endl ; 
  create_ff_con_dir() ; 
  cout << con_path_ff << endl ; 
  
  cout << "make con_prob" << endl ; 
  make_con_prob_ff() ; 
  cout << "make con_vec" << endl ; 
  make_con_vec_ff() ; 
  cout << "make sparse_rep" << endl ; 
  make_con_sparse_rep_ff() ;
  
  if(IF_SAVE_SPARSE_REP_FF) {
    cout << "save sparse_rep" << endl ; 
    save_ff_to_con_file() ; 
  }
}

void get_con_sparse_vec_ff() { 

  if(IF_SPEC || IF_RING || IF_GAUSS) {
    cout << "init theta_ff" << endl ; 
    init_theta_ff() ;
  }
  
  cout << "ff connectivity path" << endl ; 
  create_ff_con_dir() ; 
  cout << con_path_ff << endl ; 
  
  n_post_ff = new int [N_POISSON]() ; 
  read_from_file(con_path_ff, "n_post_ff", n_post_ff, N_POISSON) ; 
  
  cout << "average number of postsynaptic connections per FF neuron: " << endl ; 
  for(i=0; i<n_pop; i++) { 
    cout << "pop " << i << " " ; 
    cout << n_post_ff[i] << " " ; 
    cout << endl ;     
  }
  
  idx_post_ff = new unsigned long [N_POISSON]() ; 
  read_from_file(con_path_ff, "idx_post_ff", idx_post_ff, N_POISSON) ; 
  
  cout << "postsynaptic neurons idx: " ; 
  for(i=0; i<10; i++) 
    cout << idx_post_ff[i] << " " ; 
  cout << endl ; 
  
  for(i=0 ; i<N_POISSON; i++) 
      total_n_post_ff += n_post_ff[i] ; 
  
  cout << "total number of FF connections: " << total_n_post_ff << endl ; 
  id_post_ff = (unsigned long *) malloc( (unsigned long) total_n_post_ff * sizeof(unsigned long) ) ;
  read_from_file(con_path_ff, "id_post_ff", id_post_ff, total_n_post_ff) ;
  
  cout << "postsynaptic neurons id: " ; 
  for(int i=0; i<10; i++) 
    cout << id_post_ff[i] << " " ; 
  cout << endl ; 
}

#endif
