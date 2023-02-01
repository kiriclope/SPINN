#ifndef _UTILS_ 
#define _UTILS_

///////////////////////////////////////////////////////////////////////////
// Recurrent connectivity
///////////////////////////////////////////////////////////////////////////

void init_con_globals() { 
  cout << "initialize globals" << endl ; 
  
  IF_STRUCTURE = IF_RING || IF_SPEC || IF_LOW_RANK || IF_GAUSS ; 
  
  n_per_pop = new unsigned long [n_pop]() ;
  
  for(i=0;i<n_pop;i++)
    n_per_pop[i] = (unsigned long) ( n_frac[i] * (float) n_neurons ) ;
  
  cum_n_per_pop = new unsigned long [n_pop+1]() ;
  
  for(i=0;i<n_pop+1;i++)
    for(j=0;j<i;j++)
      cum_n_per_pop[i] += n_per_pop[j] ;
  
  which_pop = new int [n_neurons] ;
  
  for(i=0;i<n_pop;i++)
    for(j=0; j<n_neurons; j++)
      if(j>=cum_n_per_pop[i] && j<cum_n_per_pop[i+1])
  	which_pop[j] = i ;

  K_over_Na = new float [n_pop]() ; 
  
  for(i=0; i<n_pop; i++) {
    if(IF_MATO_K)
      K_over_Na[i] = (float) ( K * n_frac[i] / (float) n_per_pop[i] ) ;
    else
      K_over_Na[i] = (float) ( K / (float) n_per_pop[i] ) ;
    cout << K_over_Na[i] << " " ;
  }
  cout << endl ;
  
  if(IF_GAUSS) 
    prefactor = new float [n_pop*n_neurons]() ;
  
}

void delete_con_globals() { 
  cout << "delete globals" << endl ; 
  
  delete [] n_per_pop ; 
  delete [] cum_n_per_pop ; 
  delete [] which_pop ; 
  
  /* delete [] con_prob ;  */
  /* delete [] con_vec ;  */
  
  delete [] n_post ; 
  delete [] idx_post ;
  free(id_post) ; 
    
  if(IF_GAUSS) {
    delete [] prefactor ; 
    delete [] X ; 
  }
  
}

////////////
// Gauss
////////////

void init_X() { 
  
  X = new float [n_neurons]() ; 
  
  for(i=0; i<n_pop; i++) 
    for(j=cum_n_per_pop[i]; j<cum_n_per_pop[i+1]; j++)      
      X[j] = fmod( (float) (j-cum_n_per_pop[i]) , (float) n_per_pop[i] ) * 2.0 * M_PI / (float) n_per_pop[i] ;     
  
  cout << "X: " ; 
  for(i=0; i<10; i++) 
    cout << X[i] << " " ; 
  cout << endl ; 
  
}

float Gaussian1D(float mu, float sigma) {
  
  mu = (float) min( (float)  abs(mu), (float) (2.0*M_PI-abs(mu) ) ) ;  // mato et al.
  sigma *= DEG_TO_RAD ;
  
  if(sigma!=0.) 
    return exp(-mu*mu/2./sigma/sigma)/sqrt(2.*M_PI)/sigma ; 
  else 
    return 1. ; 
}

////////////
// Ring
////////////
void init_theta() {
  
  theta = new float [n_neurons]() ; 
  
  for(i=0; i<n_pop; i++) { 
    for(j=cum_n_per_pop[i]; j<cum_n_per_pop[i+1]; j++) 
      theta[j] = ( 2.0 * M_PI * (float) (j-cum_n_per_pop[i]) ) / (float) n_per_pop[i] ; 
  } 
  
  cout << "theta: " ; 
  for(i=0; i<10; i++) 
    cout << theta[i] << " " ; 
  cout << endl ;
  
}

////////////

void func_con_prob() {
  
  K_over_Na = new float [n_pop]() ; 
  
  for(i=0; i<n_pop; i++) {
    if(IF_MATO_K)
      K_over_Na[i] = (float) ( K * n_frac[i] / (float) n_per_pop[i] ) ;
    else
      K_over_Na[i] = (float) ( K / (float) n_per_pop[i] ) ;
    cout << K_over_Na[i] << " " ;
  }
  cout << endl ;
    
  con_prob = new float [n_neurons*n_neurons]() ; 
  
  if(IF_STRUCTURE) {
    
    if(!IF_LOW_RANK)
      init_theta() ; 
    
    if(IF_SPEC)
      kappa = KAPPA/sqrt_Ka[0] ; 
    else
      kappa = KAPPA ; 
    
    if(IF_RING || IF_SPEC) {
      
      float cos_D_theta=0.0 ; 

      if(IF_RING) cout << "ring: " ; 
      if(IF_SPEC) cout << "random network with specific connections: rank, " << RANK << ", " ;      
      
      cout << "kappa, " << KAPPA ;
      
      if(KAPPA > 0.5) {
	cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl ;
	cout << "ERROR: KAPPA TOO LARGE" << endl ; 
	cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl ;
      }            
      
      for(j=0; j<n_neurons; j++) { //pre
	pre_pop = which_pop[j] ; 
	
	kappa_K_N = K_over_Na[pre_pop] * kappa ; 
	
	for(i=0; i<n_neurons; i++) { // post
	  post_pop = which_pop[i] ; 
	  
	  con_prob[j + i * n_neurons] += K_over_Na[pre_pop] ; 
	  
	  if(IS_STRUCT_SYN[pre_pop + post_pop * n_pop]) { 
	    cos_D_theta = cos(theta[i]-theta[j]) ;	    
	    con_prob[j + i * n_neurons] += kappa_K_N * cos_D_theta ;
	    
	    /* if(RANK==2)  */
	    /*   con_prob[idx_perm[j] + idx_perm[i] * n_neurons] += kappa_1_K_N * cos_D_theta ;  */

	    if(con_prob[j + i * n_neurons]<=0 || con_prob[j + i * n_neurons]>=1) {
	      cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl ;
	      cout << "error con_prob>1 or <0"  << endl ;
	      cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl ;
	    }
	  }	  
	} 
      } 
    } //endif ring/spec 
    
    //////////
    // GAUSS
    //////////    
    if(IF_GAUSS) { 
      cout << "gaussian structure: " ;
      cout << " sigmas " ; 
      cout << SIGMA[0] <<" "<< SIGMA[1] ; 
      cout <<" "<< SIGMA[2] <<" "<< SIGMA[3] <<  endl ; 
      
      /* init_X() ;  */
      init_theta() ;
      
      for(i=0; i<n_neurons; i++) { // post, need loop over post before pre to compute prefactor
	post_pop = which_pop[i] ; 
	
	for(j=0; j<n_neurons; j++) { // pre, Jij: j (pre) to i (post)
	  pre_pop = which_pop[j] ;
	  
	  con_prob[j + i * n_neurons] = Gaussian1D(theta[i]-theta[j], SIGMA[pre_pop+post_pop*n_pop]) ; 
	  // Pij = Zb Cij with Zb = K / sum_j Cij then sum_j Pij = K 
	  prefactor[i + pre_pop*n_neurons] += con_prob[j + i * n_neurons] ;
	  // sum over presynaptic neurons j
	}

	for(j=0; j<n_neurons; j++) {
	  pre_pop = which_pop[j] ; 
	  con_prob[j + i * n_neurons] *= Ka[pre_pop] / prefactor[i+pre_pop*n_neurons] ;
	  
	  if(con_prob[j + i * n_neurons]<=0 || con_prob[j + i * n_neurons]>=1) {
	    cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl ;
	    cout << "error con_prob>1 or <0"  << endl ; 
	    cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl ;
	  }
	}
	/* prefactor[pre_pop+post_pop*n_pop] = 0.0 ; // reinit prefactor for each post  */
	 
      } // end loop post
      
    }// endif Gauss
      
    if(IF_LOW_RANK) { 
      cout << "random network with low-rank specific connections: rank, " << RANK << ", " ; 
      
      for(j=0; j<n_neurons; j++) { // pre
	pre_pop = which_pop[j] ;

	kappa_K_N = K_over_Na[pre_pop] * kappa ;
	if(RANK==2)
	  kappa_1_K_N = K_over_Na[pre_pop] * kappa_1 ; 
    
	for(i=0; i<n_neurons; i++) { // post
	  post_pop = which_pop[i] ;
	  con_prob[j + i * n_neurons] += K_over_Na[pre_pop] ;

	  if(pre_pop==0 && post_pop==0) 
	    if(IS_STRUCT_SYN[pre_pop + post_pop * n_pop]) {
	      con_prob[j + i * n_neurons] += kappa_K_N * ksi[i] * ksi[j] ;
	      if(RANK==2)
		con_prob[j + i * n_neurons] += kappa_1_K_N * ksi_1[i] * ksi_1[j] ;
	      con_prob[j + i * n_neurons] = cut_LR(con_prob[j + i * n_neurons]) ;
	    }
	  
	  if(con_prob[j + i * n_neurons]<0 || con_prob[j + i * n_neurons]>1) { 
	    cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl ; 
	    cout << "error con_prob>1 or <0"  << endl ; 
	    cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl ; 
	  }

	}

      }
      
    } // endif Low-rank 

  }// endif structure   
  else {
    cout << "random network" << endl ; 
    
    for(j=0; j<n_neurons; j++) { // Pij -> j (cols, pre) to i (rows, post) 
      pre_pop = which_pop[j] ; 
      for(i=0; i<n_neurons; i++) 
	con_prob[j + i * n_neurons] = K_over_Na[pre_pop] ; // Pij -> P[j + i * n_cols] (index = indexX * arrayWidth + indexY) ; 
    } 
  } 
  
  delete [] K_over_Na ; 
  
}

void func_con_vec() {
  
  cout << "generate connectivity vector" << endl ; 
  
  con_vec = new int [n_neurons*n_neurons]() ; 
  
  for(i=0; i<n_neurons; i++) { 
    for(j=0; j<n_neurons; j++) { // Jij -> j (cols, pre) to i (rows, post) 
      if(con_prob[j + i * n_neurons] >= unif(con_gen) ) { 
  	con_vec[j + i * n_neurons] = 1 ; 
  	total_n_post++ ; 
      } 
    } 
  } 
  
  delete [] con_prob ; 
  
  if(IF_SAVE_CON_VEC && !IF_TRIALS) { 
    cout << "saving connectivity vector to: "; 
    write_to_file(con_path, "con_vec", con_vec , n_neurons*n_neurons) ; 
  } 
} 

void func_con_sparse_rep() { 
  
  unsigned long counter = 0 ;  
  
  cout << "generate sparse representation" << endl ;
  cout << "total_n_post " << total_n_post << endl ; 
  
  idx_post = new unsigned long [n_neurons]() ; 
  id_post = (unsigned long *) malloc( (unsigned long) total_n_post * sizeof(unsigned long) ) ; 
  
  n_post = new int [n_neurons]() ;
    
  avg_n_post = new int [n_pop*n_pop]() ;
  
  /* n_pre = new int *[n_pop]() ; */
  /* for(i=0;i<n_pop;i++) // presynaptic population b */
  /*   n_pre[i] = new int [n_neurons]() ; */
  
  /* avg_n_pre = new int [n_pop*n_pop]() ;  */
  
  for(j=0; j<n_neurons; j++) { // Jij -> j (cols, pre) to i (rows, post) 
    pre_pop = which_pop[j] ; 
    
    for(i=0; i<n_neurons; i++) { 
      post_pop = which_pop[i] ;
      
      if(con_vec[j + i * n_neurons]) { // j->i = 1 with proba Kj/Nj 
	id_post[counter] = i ; 
	n_post[j]++ ; // Kb to post pop a, K/N * N = K 
	/* n_pre[pre_pop][i]++ ; // Ka from pre pop a, K/N * Nj = Kj */ 
	avg_n_post[pre_pop + post_pop*n_pop]++ ; 
	counter++ ; 
      } 
    } 
  } 
  
  delete [] con_vec ; 
  
  cout << "average number of postsynaptic connections per pop: " ;
  for(i=0; i<n_pop; i++)
    for(j=0; j<n_pop; j++)
      cout << i << j << " " << avg_n_post[j+i*n_pop] / n_per_pop[i] << " " ;
  cout << endl ;
  
  /* for(i=0; i<2*n_pop; i++) */
  /*   cout << avg_n_post[i] / n_per_pop[i%2] << " " ; */
  /* cout << endl ; */
  
  delete [] avg_n_post ;
  
  /* cout << "number of presynaptic connections per neuron: " << endl ; */

  /* for(i=0; i<n_pop; i++) { //post */
  /*   for(j=0; j<n_pop; j++) { //pre */
  /*     cout << "neuron in " << i << " n_pre in " << j << " | " ; */
  /*     for(k=n_per_pop[i]-1; k>n_per_pop[i]-10; k--)  */
  /* 	cout << n_pre[j][k] << " " ; */
  /*     cout << endl ; */
  /*   } */
  /* } */
  
  /* delete [] n_pre ; */
  
  /* cout << "number of postsynaptic connections per neuron: " << endl ;  */
  
  /* for(i=0; i<n_pop; i++) { //post  */
  /*   for(j=0; j<n_pop; j++) { //pre  */
  /*     cout << "neuron in " << j << " n_post in " << i << " " ;  */
  /*     for(k=cum_n_per_pop[j]; k<cum_n_per_pop[j]+10; k++)  */
  /* 	cout << n_post2[i][k] << " " ;  */
  /*     cout << endl ;  */
  /*   }  */
  /* }  */
   
  cout << "total number of postsynaptic connections per neuron: " << endl ;
  
  for(j=0; j<n_pop; j++) { //pre
    cout << j << " " ;
    for(k=n_per_pop[j]; k>n_per_pop[j]-10; k--)
      cout << n_post[k] << " " ; 
    cout << endl ;
  }
  
  idx_post[0] = 0 ; 
  for(i=1; i<n_neurons; i++) 
    idx_post[i] = idx_post[i-1] + n_post[i-1] ; 
  
}

void check_sparse_rep() { 
  
  cout << "checking sparse representation" << endl ; 
  con_vec = new int [n_neurons*n_neurons]() ; 
  
  for(j=0; j<n_neurons; j++) // Jij -> j (cols, pre) to i (rows, post) 
    for(i=idx_post[j]; i<idx_post[j] + n_post[j]; i++) 
      con_vec[j + id_post[i] * n_neurons] = 1 ; 
  
  cout << "saving rebuild connectivity vector to: "; 
  write_to_file(con_path, "con_vec", con_vec , n_neurons*n_neurons) ; 
  
  delete [] con_vec ;
}

void create_con_dir() { 
  
  con_path += to_string(n_pop)+"pop"; 
  
  ostringstream str_seed;  
  str_seed << SEED_CON ; 
  
  if(IF_CON_DIR)
    con_path += "/seed_" + str_seed.str() ;
  
  if(n_pop==1)
    con_path += "/N" + to_string(n_per_pop[0]/1000) ; 
  else 
    con_path += "/NE_" + to_string(n_per_pop[0]/1000) +  "_NI_" + to_string(n_per_pop[1]/1000) ; 
  
  con_path += "/K" + to_string((int)K) ;
  
  ostringstream str_kappa, str_kappa_1 ;
  str_kappa << fixed << setprecision(2) << KAPPA ; 
  str_kappa_1 << fixed << setprecision(2) << KAPPA_1 ; 
  
  ostringstream str_map_seed ; 
  str_map_seed << fixed << setprecision(0) << MAP_SEED ; 

  ostringstream str_EE, str_EI, str_IE, str_II ;  
  str_EE << fixed << setprecision(0) << SIGMA[0] ; 
  str_EI << fixed << setprecision(0) << SIGMA[1] ; 
  str_IE << fixed << setprecision(0) << SIGMA[2] ; 
  str_II << fixed << setprecision(0) << SIGMA[3] ; 

if(IF_LOW_RANK) { 
  ostringstream str_seed_ksi ; 
  str_seed_ksi << fixed << setprecision(0) << SEED_KSI ; 

  con_path += "/low_rank/rank_" + to_string(RANK) ; 
 
  if(FIX_KSI_SEED) 
    con_path += "/seed_ksi_" + str_seed_ksi.str() ; 
}
  if(IF_RING)
    con_path += "/ring/kappa_" + str_kappa.str() ; 
    
  if(IF_SPEC) { 
    con_path += "/spec/kappa_" + str_kappa.str() ; 
    if(RANK==2) {
      con_path += "_kappa_1_" + str_kappa_1.str() ; 
      if(FIX_MAP_SEED) 
	con_path += "/seed_" + str_map_seed.str() ; 
    }

  }

  if(IF_GAUSS) {
      con_path += "/gauss/EE_" + str_EE.str() ;
      con_path += "EI_" + str_EI.str() ;
      con_path += "IE_" + str_IE.str() ;
      con_path += "II_" + str_II.str() ;
  }    
    
    /* if(IF_INI_COND)  */
    /*   con_path += "/ini_cond_" + to_string( (int) INI_COND_ID ) ;  */
    
    /* if(IF_TRIALS)  */
    /*   con_path += "/trial_" + to_string( (int) TRIAL_ID ) ;  */

    make_dir(con_path) ; 
    cout << "created directory: " ; 
    cout << con_path << endl ; 

}

void save_to_con_file() {
    
  cout << "save to file" << endl ;
    
  write_to_file(con_path, "id_post", id_post, total_n_post) ; 
  write_to_file(con_path, "idx_post", idx_post, n_neurons) ; 
  write_to_file(con_path, "n_post", n_post, n_neurons) ; 
  
}

#endif
