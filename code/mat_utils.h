#ifndef __MATRIXUTILS__ 
#define __MATRIXUTILS__ 

void get_con_sparse_vec() { 

  if(IF_SPEC || IF_RING || IF_GAUSS) 
    init_theta() ;
    
  create_con_dir() ; 
  
  cout << "###############################################" << endl ; 
  cout << "getting connectivity sparse vectors from:" ; 
  cout << con_path << endl ; 
  
  n_post = new int [n_neurons]() ; 
  read_from_file(con_path, "n_post", n_post, n_neurons) ; 
  
  cout << "average number of postsynaptic connections per neuron: " << endl ; 
  for(i=0; i<n_pop; i++) { 
    cout << "pop " << i << " " ; 
    for(j=cum_n_per_pop[i]; j<cum_n_per_pop[i]+10; j++) 
      cout << n_post[j] << " " ; 
    cout << endl ;     
  }
  
  idx_post = new unsigned long [n_neurons]() ; 
  read_from_file(con_path, "idx_post", idx_post, n_neurons) ; 
  
  cout << "postsynaptic neurons idx: " ; 
  for(i=0; i<10; i++) 
    cout << idx_post[i] << " " ; 
  cout << endl ; 
  
  for(i=0 ; i<n_neurons; i++) 
      total_n_post += n_post[i] ; 
  
  cout << "total number of connections: " << total_n_post << endl ; 
  id_post = (unsigned long *) malloc( (unsigned long) total_n_post * sizeof(unsigned long) ) ;
  /* id_post = new unsigned long [total_n_post]() ; // for some reason this is not working ... */ 
  read_from_file(con_path, "id_post", id_post, total_n_post) ;
  
  cout << "postsynaptic neurons id: " ; 
  for(int i=0; i<10; i++) 
    cout << id_post[i] << " " ; 
  cout << endl ;

  if(IF_LOW_RANK) 
    get_ksi() ; 
  
} 


void gen_con_sparse_vec() { 
  cout << "###############################################" << endl ;

  cout << "connectivity path" << endl ; 
  create_con_dir() ;
  
  cout << "make con_prob" << endl ;  
  func_con_prob() ;
  cout << "make con_vec" << endl ; 
  func_con_vec() ;
  cout << "make sparse_rep" << endl ; 
  func_con_sparse_rep() ;
  
  if(IF_CHECK_SPARSE_REP) 
    check_sparse_rep() ; 

  if(IF_SAVE_SPARSE_REP) 
    save_to_con_file() ; 
  
}

#endif
