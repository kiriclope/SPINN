#ifndef __HYST__ 
#define __HYST__ 

void hyst_init() { 
  
  cout << "hyst_min " << HYST_X_MIN << " hyst_max " << HYST_X_MAX << " hyst_dX " << HYST_DX << endl ; 
  
  duration = (HYST_X_MAX - HYST_X_MIN + 2.0 * HYST_DX ) / HYST_DX * ( TIME_WINDOW + TIME_STEADY ) ; 
  time_rec = duration ; 
  time_rec_spikes = 2000.0 ; 
  
  if(HYST_J_EE==1) { 
    if(IF_LIF) 
      J_scaled[0] = GAIN / sqrt_Ka[0] / TAU_SYN[0] * (Vth-Vr[0]) * HYST_X_MIN ; 
    if(IF_BIN) 
      J_scaled[0] = GAIN / sqrt_Ka[0] * HYST_X_MIN ; 
    
    J[0] = HYST_X_MIN ; 
  } 
  
  if(HYST_J_EE==-1) { 
    if(IF_LIF) 
      J_scaled[0] = GAIN / sqrt_Ka[0] / TAU_SYN[0] * (Vth-Vr[0]) * HYST_X_MAX ; 
    if(IF_BIN) 
      J_scaled[0] = GAIN / sqrt_Ka[0] * HYST_X_MAX ; 
    
    J[0] = HYST_X_MAX ; 
  }
  
  if(HYST_M0==1) { 
    
    for(i=0;i<n_pop;i++) { 
      m0 = HYST_X_MIN ; 
      if(IF_LIF) 
	ext_inputs_scaled[i] = ext_inputs[i] * sqrt_Ka[0] * (Vth-Vr[0]) * m0 ; 
      if(IF_BIN) 
	ext_inputs_scaled[i] = ext_inputs[i] * sqrt_Ka[0] * m0 ; 
    } 
    
    for(i=0;i<n_neurons;i++) 
      ff_inputs[i] = ext_inputs_scaled[which_pop[i]] ;    
  }
  
  if(HYST_M0==-1) {
    
    for(i=0;i<n_pop;i++) { 
      m0 = HYST_X_MAX ; 
      if(IF_LIF) 
	ext_inputs_scaled[i] = ext_inputs[i] * sqrt_Ka[0] * (Vth-Vr[0]) * m0 ; 
      if(IF_BIN) 
	ext_inputs_scaled[i] = ext_inputs[i] * sqrt_Ka[0] * m0 ; 
    } 
    
    for(i=0;i<n_neurons;i++) 
      ff_inputs[i] = ext_inputs_scaled[which_pop[i]] ; 
  }

  cout << "m0 " << m0 << " Jee " << J[0] << endl ; 
} 

void hyst_update() { 

  if(HYST_J_EE!=0) {
    J[0] += HYST_DX * HYST_J_EE ;
    
    if(IF_LIF) 
      J_scaled[0] = GAIN * J[0] / sqrt_Ka[0] / TAU_SYN[0] * (Vth-Vr[0]) ;
    
    if(IF_BIN) 
      J_scaled[0] = GAIN * J[0] / sqrt_Ka[0] ; 
    
    if(J[0] > HYST_X_MAX + HYST_DX || J[0] < HYST_X_MIN - HYST_DX ) 
      exit(0) ; 
  }
  
  if(HYST_M0!=0) {
    m0 += HYST_DX * HYST_M0 ; 
    
    if(IF_LIF) 
      for(i=0;i<n_pop;i++) 
	ext_inputs_scaled[i] = ext_inputs[i] * sqrt_Ka[0] * (Vth-Vr[0]) * m0 ;
    
    if(IF_BIN) 
      for(i=0;i<n_pop;i++) 
	ext_inputs_scaled[i] = ext_inputs[i] * sqrt_Ka[0] * m0 ; 
    
    for(i=0;i<n_neurons;i++) 
      ff_inputs[i] = ext_inputs_scaled[which_pop[i]] ; 
      
    if(m0 > HYST_X_MAX + HYST_DX || m0 < HYST_X_MIN - HYST_DX ) 
      exit(0) ; 
  } 
  
  time_steady += TIME_STEADY + TIME_WINDOW ; 

  /* cout << "m0 " << m0 << " Jee " << J[0] << " Jee/sqrt_K "<< J_scaled[0] ;
     cout << " Ie "<< ext_inputs[0] << " sqrt_K Ie " << ext_inputs_scaled[0] << endl ; */ 
  
}     

#endif 
