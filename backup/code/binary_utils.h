#ifndef __BINUTILS__ 
#define __BINUTILS__ 

#include <vector>
#include <functional> // std::plus 
#include <numeric> //std::accumulate

vector<int> spin, cum_spin, n_updates ; 
int n_updates_tot=0 ; 
const double THETA[2] = {1.0, 1.0} ; 
const double TAU[2] = {1.0, 10.0} ; 

double dt ; 
double proba[2] ;

uniform_real_distribution<double> unif_one(0.0, 1.0) ; 
uniform_int_distribution<int> unif_pop ; 

unsigned long sum(vector<int> v,  int a, int b) {
  unsigned long sum = 0 ; 
  for(unsigned long k=a; k<b; k++) 
    sum += v[k] ; 
  return sum ; 
}

void init_bin_globals() {
  t_spike = new double [n_neurons]() ; // time of spike emission
  
  spin.resize(n_neurons) ; 
  cum_spin.resize(n_neurons) ; 
  n_updates.resize(n_neurons) ; 
  
  for(i=0;i<n_pop;i++) { 
    ext_inputs_scaled[i] = GAIN * ext_inputs[i] * sqrt_Ka[0] * m0 ; 
    
    for(j=0;j<n_pop;j++) 
      /* J_scaled[j+i*n_pop] = GAIN * J[j+i*n_pop] / sqrt_Ka[j] * sqrt_Ka[0] / sqrt_Ka[j] ;  */
      J_scaled[j+i*n_pop] = GAIN * J[j+i*n_pop] / sqrt_Ka[j] ; 
  }
  
  for(i=0; i<n_neurons; i++) { 
    ff_inputs[i] = ext_inputs_scaled[which_pop[i]] ; 
    net_inputs[i] = ext_inputs_scaled[which_pop[i]] ; 
  }
  
  if(n_pop==1) 
    dt =  TAU[0] / (double) n_per_pop[0] ; 
  else 
    dt = (TAU[0] * TAU[1]) / ( (double) n_per_pop[0] * TAU[1] + (double) n_per_pop[1] * TAU[0]) ; 
  
  cout << "dt " << dt << " update proba: " ; 
  for(i=0; i<n_pop; i++) { 
    proba[i] = dt * (double) n_per_pop[i] / TAU[i] ; 
    cout << proba[i] << " " ; 
  }
  cout << endl ; 
  
  duration = DURATION + TIME_STEADY + TIME_WINDOW ; 
  time_rec = TIME_REC + TIME_STEADY + TIME_WINDOW ; 
  time_rec_spikes = TIME_REC_SPIKES + TIME_STEADY ; 
  time_steady = TIME_STEADY ;   
}

void delete_bin_globals() { 
  spin.clear() ; 
  cum_spin.clear() ; 
  n_updates.clear() ; 
  delete [] t_spike ; 
}
  
void pick_pre_pop_neuron() { 

  if(n_pop==2)
    if(unif_one(rand_gen)<proba[0]) 
      pre_pop = 0 ; 
    else 
      pre_pop = 1 ; 
  else
    pre_pop = 0 ; 
  
  unif_pop.param(uniform_int_distribution<int>::param_type(0, n_per_pop[pre_pop]-1) ) ; 
  
  i_neuron = unif_pop(rand_gen) ; 
  i_neuron += cum_n_per_pop[pre_pop] ; 
}

void update_post_inputs() { 
  // only updating inputs of postsynaptic targets,  
  // if Inet>=1: Ipost+=J iif spin 0->1 elif Inet<1: Ipost-=J iif spin 1->0 (using eps={1,-1}) 
  for(j=idx_post[i_neuron]; j<idx_post[i_neuron] + n_post[i_neuron]; j++) {  
    post_pop = which_pop[id_post[j]] ;
    
    if(IF_STP && stp_synapse[pre_pop + post_pop * n_pop]) 
      inputs[pre_pop][id_post[j]] += eps[spin[i_neuron]] * A_u_x_stp[i_neuron] * J_scaled[pre_pop + post_pop * n_pop] ; 
    else 
      inputs[pre_pop][id_post[j]] += eps[spin[i_neuron]] * J_scaled[pre_pop + post_pop * n_pop] ; 
  } 
} 

void update_spins_and_post_inputs() { 
  
  net_inputs[i_neuron] = ff_inputs[i_neuron] ; 
  for(i=0; i<n_pop; i++) 
    net_inputs[i_neuron] += inputs[i][i_neuron] ; 
  
  if(net_inputs[i_neuron]>=THETA[pre_pop]) { // if spike

    n_updates[i_neuron]++ ; 
    n_updates_tot++ ; 
    
    if(spin[i_neuron]==0) { 
      
      ISI = t_time - t_spike[i_neuron] ;
      t_spike[i_neuron] = t_time ;
      
      if(IF_STP) 
	update_stp_variables_lif() ; 
      
      update_post_inputs() ; 
      spin[i_neuron] = 1 ; 
      
      if(t_time>time_steady && t_time<time_rec_spikes) // save spike times 
	file_spike_times << fixed << setprecision(1) << (float) (i_neuron) 
			 << " " << (float) (t_time-TIME_STEADY) << endl ;
      
    } 
  } 
  else { 
    if(spin[i_neuron]==1) { 
      update_post_inputs() ; 
      spin[i_neuron] = 0 ; 
    } 
  } 
} 

void initial_conditions_bin() { 
  // THIS IS MESSING UP WITH THE STEADY STATE DO NOT USE !!!!!!!!
  normal_distribution<double> normal(0.0, 1.0) ; 
  uniform_real_distribution<double> unif(0.0, 1.0) ; 
  
  if(IF_STP) {
    
    for(i=0;i<n_neurons;i++) {
      x_stp[i] = unif(rand_gen) ;
      u_stp[i] = unif(rand_gen) ;
    }
    
  }
  
  /* double mean, sigma ;  */
  
  /* for(i=0;i<n_pop;i++) {  */
  /*   mean = unif(rand_gen) ;  */
  /*   sigma = unif(rand_gen) ; // divide by 4 so that the distribution is between 0 and 1 at best  */
    
  /*   for(i_neuron=0;i_neuron<n_neurons;i_neuron++)  */
  /*     if(i==0)  */
  /* 	inputs[i][i_neuron] = max( mean + sqrt(sigma) * normal(rand_gen), 0.0 ) ;  */
  /*     else  */
  /* 	inputs[i][i_neuron] = min(-mean + sqrt(sigma) * normal(rand_gen), 0.0 ) ;  */
  /* }  */
  
  /* for(i_neuron=0; i_neuron<n_neurons; i_neuron++) {  */
  /*   net_inputs[i_neuron] = ff_inputs[i_neuron] ;  */
  /*   for(i=0; i<n_pop; i++)  */
  /*     net_inputs[i_neuron] += inputs[i][i_neuron] ;  */
  /* } */
  
  // this is not working 
  /* for(i_neuron=0; i_neuron<n_neurons; i_neuron++) { */ 
  
  /*   if(net_inputs[i_neuron]>=THETA[pre_pop]) { // if spike */ 
  
  /*     ISI = t_time - t_spike[i_neuron] ; */ 
  /*     t_spike[i_neuron] = t_time ; */ 
  
  /*     if(IF_STP) */ 
  /* 	update_stp_variables_lif(ISI) ; */ 
  
  /*     update_post_inputs() ;  */ 
  /*     spin[i_neuron] = 1 ;  */ 
  
  /*   } */ 
  /* } */ 
  
  double ini_proba[2] ; 

  if(n_pop==1)
    ini_proba[0] = 0.5 ;
  else
    while(ini_proba[0]>=0.75*ini_proba[1]) {
      ini_proba[0] = unif(rand_gen) ;
      ini_proba[1] = unif(rand_gen) ;
    }
  
  cout << "ini_proba " << ini_proba[0] << " " << ini_proba[1] << endl ; 
  
  for(i_neuron=0; i_neuron<n_neurons; i_neuron++) {
    
    pre_pop = which_pop[i_neuron] ;
    
    if(ini_proba[pre_pop]>=unif(rand_gen)) { 
      
      ISI = t_time - t_spike[i_neuron] ; 
      t_spike[i_neuron] = t_time ; 
      
      if(IF_STP) 
  	update_stp_variables_lif() ; 
      
      update_post_inputs() ; 
      spin[i_neuron] = 1 ; 
    } 
    
  } 
  
  /* cout << "initial activities: " ;  */
  /*   for(i=0;i<n_pop; i++)  */
  /*     cout << " " << (double) sum(spin, cum_n_per_pop[i], cum_n_per_pop[i+1]) / n_per_pop[i] ;  */
  /*   cout << "\n" ;      */
  
  /* double phi_ini = unif(rand_gen) * 2.0 * M_PI ; */
  /* double kappa_ini = unif(rand_gen) ; */
  /* cout << "kappa_ini " << kappa_ini << " phi ini " << phi_ini << endl ; */
  
  /* for(i=0;i<n_per_pop[0]; i++) { */
  /*   /\* ff_inputs[i] = ext_inputs_scaled[which_pop[i]] * ( 1.0 + 1.0 / sqrt_Ka[which_pop[i]] * phi_ini ) ;  *\/ */
    
  /*   ff_inputs[i] = ext_inputs_scaled[0] * ( 1.0 +  kappa_ini / sqrt_Ka[0]  */
  /*   						       * cos( 2.0 * (double) i * M_PI / (double) n_per_pop[0] + phi_ini) ) ;  */
    
  /*   /\* ff_inputs[i] = ext_inputs_scaled[0] * ( 1.0 + unif(rand_gen) / sqrt_Ka[0] ) ;  *\/  */
  /* } */
  
  /* for(t_time=0.0; t_time<TIME_INI; t_time+=dt) {  */
  /*   pick_pre_pop_neuron() ;  */
  /*   update_spins_and_post_inputs() ;  */
  /* }// endfor time loop  */
  
  /* for(i=0;i<n_per_pop[0]; i++)  */
  /*   ff_inputs[i] = ext_inputs_scaled[0] ; */
  
  /* for(i_neuron=0; i_neuron<n_neurons; i_neuron++)  */
  /*   t_spike[i_neuron] -= TIME_INI ;  */
  
} 

void print_activities() { 
  
  if( t_window<TIME_WINDOW-dt ) { 
    cout << fixed << setprecision(2) << round( percentage*1000.0 )/10.0 << "% " ; 
    cout << "time " << t_time << " ms " << "<n_updates> " << (double) n_updates_tot /(double) n_neurons << " " ;
    cout << "\r" ; 
  }
  
  if(t_window>=TIME_WINDOW-dt) { 
    
    cout << fixed << setprecision(2) << round( percentage*1000.0 )/10.0 << "% " ; 
    
    cout << "time " << t_time << " ms " << "<n_updates> " << (double) n_updates_tot /(double) n_neurons << " " ;
    if(HYST_J_EE!=0) 
      cout << " J_EE " << J[0] ; 
    if(HYST_M0!=0 || IF_LOOP_M0) 
      cout << "| m0 " << m0 ; 
    
    for(i=0;i<n_pop;i++)
      mean_rates[i] = (double) sum(cum_spin, cum_n_per_pop[i], cum_n_per_pop[i+1]) / n_per_pop[i] * dt / TIME_WINDOW ; 
    
    cout << " | activities:" ; 
    for(i=0;i<n_pop; i++) 
      cout << " " <<  mean_rates[i] ; 

    for(i=0; i<n_neurons;i++)
      filter_rates[i] = (double) cum_spin[i] * dt / TIME_WINDOW ;

    if(IF_SPEC) {
      get_m1_phase() ;
      cout << " m1 " ; 
      for(i=0;i<n_pop;i++)
	cout << m1[i] << " " ;
      
      cout << "phase " ;       
      for(i=0;i<n_pop;i++)
	cout << phase[i] << " " ;
    }
    
    if(IF_LOW_RANK) { 
      cout << " | overlaps:" ; 
      for(i=0;i<n_pop;i++) 
	cout << " " << overlaps[i] * 1000. / TIME_WINDOW * IS_STRUCT_SYN[i] ; 
    } 
    
    cout << "\r" ;
  }
  
}

void save_bin_file() { 
    
  if(IF_HYSTERESIS==0)
    file_mean_rates << t_time - time_steady ;
  else {
    if(HYST_J_EE!=0) 
      file_mean_rates << J[0] ; 
    if(HYST_M0!=0 || IF_LOOP_M0) 
      file_mean_rates << m0 ; 
  } 
  
  for(i=0;i<n_pop; i++) 
    file_mean_rates << " " << mean_rates[i] ; 
  file_mean_rates << endl ; 
  
  if(IF_HYSTERESIS==0) 
    file_filter_rates << t_time - time_steady ;
  else {
    if(HYST_J_EE!=0)
      file_filter_rates << J[0] ;
    if(HYST_M0!=0 || IF_LOOP_M0) 
      file_filter_rates << m0 ; 
  }
  
  for(i=0; i<n_neurons; i++) 
    file_filter_rates << " " << filter_rates[i] ; 
  file_filter_rates << endl ; 

} 

void run_sim_bin() { 

  open_files() ;
  
  init_bin_globals() ;   
  open_lif_files() ; 

  if(IF_STP) {
    open_stp_files() ;
    init_stp_globals() ;
  }
  
  initial_conditions_bin() ;
  
  if(IF_HYSTERESIS) 
    hyst_init() ; 
  
  print_sim_info() ; 

  for(t_time=0.0; t_time<duration; t_time+=dt) { 

    percentage = t_time/duration ; 
    
    if(IF_STIM) 
      tasks_inputs() ; 
    
    pick_pre_pop_neuron() ; 
    /* cout << "pre_pop "<< pre_pop << " i_neuron " << i_neuron << " net_inputs " << net_inputs[i_neuron] << endl ; */ 
    update_spins_and_post_inputs() ; 
    
    print_activities() ; 
    
    if(t_time>=time_steady) { 
      t_window += dt ; 
      transform(cum_spin.begin(), cum_spin.end(), spin.begin(), cum_spin.begin(), plus<int>()) ; 
      
      if(t_window>=TIME_WINDOW) {
	
      	save_bin_file() ; 
	fill(cum_spin.begin(), cum_spin.end(), 0) ; 
	
	if(IF_HYSTERESIS) 
	  hyst_update() ; 
	
	if(IF_STP) 
	  save_xy_to_file() ; 
	
	t_window = 0 ; 
      } 
    } 
    
  }// endfor time loop 
  
  delete_bin_globals() ; 
  delete_globals() ; 
  
  close_files() ; 
  close_lif_files() ; 
  
  if(IF_STP) { 
    delete_stp_globals() ; 
    close_stp_files() ; 
  }
  
}

#endif
