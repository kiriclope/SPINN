#ifndef __TASKSUTILS__ 
#define __TASKSUTILS__ 

void track_input() {
  float dum1=0, dum2=0, dum3=0, dum4=0 ; 
    
  if(t_time-TIME_STEADY > 2000.0 && dum1==0) {
    dum1=1 ;
    for(i=0;i<n_per_pop[0]; i++) 
      ff_inputs[i] = ext_inputs_scaled[which_pop[i]] * ( 1.0 +  KAPPA_EXT / sqrt(K) * cos( 2.0 * (float) i * M_PI / (float) n_per_pop[which_pop[i]] ) ) ; 
  }
  
  if(t_time-TIME_STEADY > 6000.0 && dum2==0) { 
    dum2=1 ; 
    for(i=0;i<n_per_pop[0]; i++) 
      ff_inputs[i] = ext_inputs_scaled[which_pop[i]] * ( 1.0 +  KAPPA_EXT / sqrt(K) * cos( 2.0 * (float) i * M_PI / (float) n_per_pop[which_pop[i]] + M_PI/2.0) ) ; 
  }

  if(t_time-TIME_STEADY > 10000.0 && dum3==0) { 
    dum3=1 ; 
    for(i=0;i<n_per_pop[0]; i++) 
      ff_inputs[i] = ext_inputs_scaled[which_pop[i]] * ( 1.0 +  KAPPA_EXT / sqrt(K) * cos( 2.0 * (float) i * M_PI / (float) n_per_pop[which_pop[i]] + M_PI) ) ; 
  }
    
  if(t_time-TIME_STEADY > 14000.0 && dum4==0) {
    dum4=1 ;
    for(i=0;i<n_per_pop[0]; i++) 
      ff_inputs[i] = ext_inputs_scaled[which_pop[i]] * ( 1.0 +  KAPPA_EXT / sqrt(K) * cos( 2.0 * (float) i * M_PI / (float) n_per_pop[which_pop[i]] + 3*M_PI/2.0 ) ) ; 
  }
  
}

void christos_task() {
  
  // CUE 
  if(t_time-TIME_STEADY >= T_CUE_ON && t_time-TIME_STEADY < T_CUE_OFF  && !SWITCH_ON) {
    
    for(i=0;i<n_neurons; i++)
      bg_inputs[i] += sqrt_Ka[0] * A_CUE[which_pop[i]] * ( 1.0 + EPS_CUE[which_pop[i]]
							   * cos( theta[i] - 2.0 * PHI_CUE * M_PI) ) ; 
    
    SWITCH_ON = 1 ; 
  }
  
  if(t_time-TIME_STEADY >= T_CUE_OFF && SWITCH_ON) { 
    
    for(i=0;i<n_neurons; i++) 
      bg_inputs[i] -= sqrt_Ka[0] * A_CUE[which_pop[i]] * ( 1.0 + EPS_CUE[which_pop[i]]
							   * cos( theta[i] - 2.0 * PHI_CUE * M_PI) ) ; 
    SWITCH_ON = 0 ; 
  } 
  
  // Second cue 
  if(t_time-TIME_STEADY >= T_ERASE_ON && t_time-TIME_STEADY < T_ERASE_OFF  && !SWITCH_OFF) {
    if(IF_SECOND)
      for(i=0;i<n_neurons; i++)
	bg_inputs[i] += sqrt_Ka[0] * A_ERASE[which_pop[i]] * ( 1.0 + EPS_ERASE[which_pop[i]]
							       * cos( theta[i] - 2.0 * PHI_ERASE * M_PI) ) ;
    else
      for(i=0;i<n_neurons; i++)
	bg_inputs[i] += A_ERASE[which_pop[i]] * ( 1.0 + EPS_ERASE[which_pop[i]]
						  * cos( theta[i] - 2.0 * PHI_ERASE * M_PI) ) ;
    
    SWITCH_OFF = 1 ;
  }
  
  if(t_time-TIME_STEADY >= T_ERASE_OFF && SWITCH_OFF) {
    if(IF_SECOND)
      for(i=0;i<n_neurons; i++)
	bg_inputs[i] -= sqrt_Ka[0] * A_ERASE[which_pop[i]] * ( 1.0 + EPS_ERASE[which_pop[i]]
							       * cos( theta[i] - 2.0 * PHI_ERASE * M_PI) ) ; 
    else
      for(i=0;i<n_neurons; i++)
	bg_inputs[i] -= A_ERASE[which_pop[i]] * ( 1.0 + EPS_ERASE[which_pop[i]]
						  * cos( theta[i] - 2.0 * PHI_ERASE * M_PI) ) ; 
    
    SWITCH_OFF = 0 ;
  }
}

void DPA_task() {
  
  // Sample stimulus ON 
  if(t_time-TIME_STEADY >= T_SAMPLE_ON && t_time-TIME_STEADY < T_SAMPLE_OFF  && !SWITCH_OFF) { 
    cout << "Sample ON " << endl ; 
    kappa_ext = KAPPA_EXT / sqrt_Ka[0] * (1.0 + white_noise(rand_gen) / 2.0 ) ; 
    
    if(IF_LOW_RANK)
      if(SAMPLE)
	for(i=0;i<n_per_pop[0]; i++) 
	  bg_inputs[i] = ext_inputs_scaled[0] * ( 1.0 + kappa_ext * sample_B[i] ) ; 
      else
	for(i=0;i<n_per_pop[0]; i++) 
	  bg_inputs[i] = ext_inputs_scaled[0] * ( 1.0 + kappa_ext * sample_A[i] ) ; 
	
    SWITCH_OFF = 1 ; 
  }
  
  // Sample stimulus OFF
  if(t_time-TIME_STEADY >= T_SAMPLE_OFF && SWITCH_OFF==1) { 
    cout << "Sample OFF " << endl ;
    for(i=0;i<n_neurons; i++) 
      bg_inputs[i] = ext_inputs_scaled[which_pop[i]] ;
    
    SWITCH_OFF = 2 ;
  }
  
  // Test stimulus ON 
  if(t_time-TIME_STEADY >= T_TEST_ON && t_time-TIME_STEADY < T_TEST_OFF && SWITCH_OFF==2) {
    kappa_test = KAPPA_TEST / sqrt_Ka[0] ; 
    cout << "Test ON " << endl ; 
    for(i=0;i<n_per_pop[0]; i++)      
      bg_inputs[i] = ext_inputs_scaled[0] * ( 1.0 +  kappa_test * ksi_1[i] ) ; 
    
    /* for(i=n_per_pop[0];i<n_per_pop[1]; i++) */
    /*   ff_inputs[i] = ext_inputs_scaled[1] * ( 1.0 +  kappa_ext * unif(rand_gen) ) ;  */

    SWITCH_OFF = 3 ; 
  } 
  
  // Test stimulus OFF 
  if(t_time-TIME_STEADY >= T_TEST_OFF && SWITCH_OFF==3) {
    cout << "Test OFF " << endl ;
    for(i=0;i<n_per_pop[0]; i++) 
      bg_inputs[i] = ext_inputs_scaled[0] ;
    
    SWITCH_OFF = 4 ; 
  }
}

void DRT_task() {
  
  // distractor ON 
  if(t_time-TIME_STEADY >= T_DIST_ON && t_time-TIME_STEADY < T_DIST_OFF  && !SWITCH_ON) { 
    kappa_dist = KAPPA_DIST / sqrt_Ka[0] * (1.0 + white_noise(rand_gen) / 2.0 ) ; 

    cout << "DIST ON " << endl ;
    
    if(IF_LOW_RANK) 
      for(i=0;i<n_per_pop[0]; i++) 
	bg_inputs[i] = ext_inputs_scaled[0] * ( 1.0 + EPS[DISTRACTOR] * kappa_dist * ksi_1[i] ) ; 
    
    if(IF_SPEC) { 
      for(i=0;i<n_per_pop[0]; i++) { 
	if(RANK==1) 
	  bg_inputs[i] = ext_inputs_scaled[0] * ( 1.0 + kappa_dist * unif(rand_gen) ) ; 
	if(RANK==2) 
	  bg_inputs[idx_perm[i]] = ext_inputs_scaled[0] * ( 1.0 + kappa_dist * cos( theta[i] + PHI_DIST) ) ; 
      }
    }
    
    SWITCH_ON = 1 ; 
  }

  // distractor OFF
  if(t_time-TIME_STEADY >= T_DIST_OFF && SWITCH_ON==1) { 
    cout << "DIST OFF " << endl ;
    for(i=0;i<n_per_pop[0]; i++) 
      bg_inputs[i] = ext_inputs_scaled[0] ; 
    SWITCH_ON = 2 ; 
  }
  
  // Cue/Reward ON
  if(t_time-TIME_STEADY >= T_RWD_ON && t_time-TIME_STEADY < T_RWD_OFF  && SWITCH_ON==2) { 
    cout << "RWD ON " << endl ;
    kappa_cue = KAPPA_CUE / sqrt_Ka[0] ; 
    
      if(IF_LOW_RANK) 
	for(i=0;i<n_per_pop[0]; i++) 
	  bg_inputs[i] = ext_inputs_scaled[0] * ( 1.0 +  kappa_cue * ksi_1[i] ) ; 
      
      if(IF_SPEC) { 
	for(i=0;i<n_per_pop[0]; i++) { 
	if(RANK==1) 
	  bg_inputs[i] = ext_inputs_scaled[0] * ( 1.0 + kappa_cue * unif(rand_gen) ) ; 
	if(RANK==2) 
	  bg_inputs[idx_perm[i]] = ext_inputs_scaled[0] * ( 1.0 + kappa_cue * cos( theta[i] + PHI_DIST) ) ; 
	}
      }
      
      SWITCH_ON = 3 ;     
  }
  
  // Rwd OFF
  if(t_time-TIME_STEADY >= T_RWD_OFF && SWITCH_ON==3) { 
    cout << "RWD OFF " << endl ;
    for(i=0;i<n_per_pop[0]; i++) 
      bg_inputs[i] = ext_inputs_scaled[0] ; 
    SWITCH_ON = 4 ; 
  }
  
} 

void christos_task_ff() {
  
  if(t_time-TIME_STEADY >= T_CUE_ON && t_time-TIME_STEADY < T_CUE_OFF  && !SWITCH_ON) {
    /* m0 = ( 1.0 + A_CUE[0] ) * M0 ; */
    /* m0 = 5. * M0 ;  */
    
    for(i=0;i<N_POISSON; i++) 
      /* J_task[i] = (1.0 + A_CUE[0]) */
      /* 	* ( 1.0 + EPS_CUE[0] * cos( theta_ff[i] - 2.0 * PHI_CUE * M_PI) ) ; */
      /* J_task[i] = A_CUE[0] * ( 1.0 + EPS_CUE[0] * cos( theta_ff[i] - 2.0 * PHI_CUE * M_PI) ) / 100.0 ; */
      J_task[i] = sqrt_Ka[0] * A_CUE[0] * ( 1.0 + EPS_CUE[0] * cos( theta[i] - 2.0 * PHI_CUE * M_PI) ) ; 
    
    SWITCH_ON = 1 ; 
  }
  
  if(t_time-TIME_STEADY >= T_CUE_OFF  && SWITCH_ON) {
    /* m0 = M0 ; */
    for(i=0;i<N_POISSON; i++) 
      J_task[i] = 0.0 ; 
    SWITCH_ON = 0 ; 
  } 


  if(t_time-TIME_STEADY >= T_ERASE_ON && t_time-TIME_STEADY < T_ERASE_OFF  && !SWITCH_OFF) {
    
    for(i=0;i<N_POISSON; i++) 
      /* J_task[i] = (1.0 + A_CUE[0]) */
      /* 	* ( 1.0 + EPS_CUE[0] * cos( theta_ff[i] - 2.0 * PHI_CUE * M_PI) ) ; */
      /* J_task[i] = A_CUE[0] * ( 1.0 + EPS_CUE[0] * cos( theta_ff[i] - 2.0 * PHI_CUE * M_PI) ) / 100.0 ; */
      J_task[i] = sqrt_Ka[0] * A_ERASE[0] * ( 1.0 + EPS_ERASE[0] * cos( theta[i] - 2.0 * PHI_ERASE * M_PI) ) ; 
    
    SWITCH_OFF = 1 ; 
  }
  
  if(t_time-TIME_STEADY >= T_ERASE_OFF  && SWITCH_OFF) {
    for(i=0;i<N_POISSON; i++) 
      J_task[i] = 0.0 ; 
    SWITCH_OFF = 0 ; 
  } 
  
}

void step_input() {
  
  // Sample stimulus ON 
  if(t_time-TIME_STEADY >= T_STEP_ON && t_time-TIME_STEADY < T_STEP_OFF  && !SWITCH_ON) {
    kappa_ext = KAPPA_EXT / sqrt_Ka[0] * (1.0 + white_noise(rand_gen) ) ; 
    
    for(i=0;i<n_per_pop[0]; i++) 
      bg_inputs[i] = ext_inputs_scaled[0] * ( 1.0 +  kappa_ext * sample_A[i] ) ; 
    SWITCH_ON = 1 ; 
  }
  
  // Sample stimulus OFF 
  if(t_time-TIME_STEADY >= T_STEP_OFF && SWITCH_ON) { 
    for(i=0;i<n_per_pop[0]; i++) 
      bg_inputs[i] =  ext_inputs_scaled[0] ; 
    SWITCH_ON = 0 ; 
  } 
} 

void tasks_inputs() { 
  if(IF_STEP)
    step_input() ;
  
  if(IF_CHRISTOS)
    if(IF_FF_LAYER)
      christos_task_ff() ;
    else
      christos_task() ; 
  
  if(IF_DPA || IF_DUAL) 
    DPA_task() ; 
  
  if(IF_DUAL || IF_DRT) 
    DRT_task() ; 
}

#endif 
