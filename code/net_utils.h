#ifndef __NETUTILS__
#define __NETUTILS__

string str_mean_rates , str_filter_rates, str_inputs, str_overlaps, str_ff_inputs ; 
ofstream file_mean_rates, file_filter_rates, file_inputs, file_overlaps, file_ff_inputs ;

int read_params() {
  Config cfg;

  // Read the file. If there is an error, report it and exit.
  try {
    cfg.readFile("params.cfg");
  }
  catch(const FileIOException &fioex)
    {
      std::cerr << "I/O error while reading file." << std::endl;
      return(EXIT_FAILURE);
    }
  catch(const ParseException &pex)
    {
      std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
		<< " - " << pex.getError() << std::endl;
      return(EXIT_FAILURE);
    }

  const Setting& root = cfg.getRoot() ; 
  
  INI_COND_ID = root["INI_COND_ID"] ; 
  TRIAL_ID = root["TRIAL_ID"] ; 

  cout << "1" << endl ;
  
  VrE = (float) root["VrE"] ;
  VlE = (float) root["VlE"] ;
  
  M0 = (float) root["M0"] ; 
  NU_FF = (float) root["NU_FF"] ;
  
  // christos
  IF_CHRISTOS = (int) root["IF_CHRISTOS"] ; 
  PHI_CUE = (float) root["PHI_CUE"] ;  
  PHI_ERASE = (float) root["PHI_ERASE"] ;
  
  // dual task 
  IF_DPA = (int) root["IF_DPA"] ;  
  IF_DUAL = (int) root["IF_DUAL"] ; 
  
  SAMPLE = (int) root["SAMPLE"] ;
  DISTRACTOR = (int) root["DISTRACTOR"] ;
  
  KAPPA = (float) root["KAPPA"] ;
  KAPPA_1 = (float) root["KAPPA_1"] ;
  
  return 0 ;
}

void print_sim_info() {
  cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% " << endl ; 

  if(IF_LIF)
    cout << "LIF model " ;
  if(IF_BIN) 
    cout << "Binary model " ; 
  if(IF_RATE)
    cout << "Rate model " ; 
      
  if(IF_RK2) 
    cout << "with RK2 scheme and interpolation " ; 
  else
    cout << "with EULER integration " ; 
      
  if(IF_STP) 
    cout << "and STP " ; 
  
  cout << endl ; 
  cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% " << endl ;

  cout << "path: " << path << endl;
  
  cout << "Main loop :" ; 
  cout << " duration " << DURATION ; 
  cout << " ms | DT " << DT ; 
  cout << " ms | TIME_STEADY " << TIME_STEADY ; 
  cout << " ms | TIME_WINDOW " << TIME_WINDOW ; 
  cout << " ms | TIME_REC " << TIME_REC ; 
  cout << " ms | TIME_REC_SPIKES " << TIME_REC_SPIKES ; 
  cout << endl ;
  
  cout << "mean field rates: " ; 
  if(IF_LIF) { 
    for(i=0; i<n_pop; i++) 
      cout << mf_rates[i] << " " ; 
    cout << endl ;
  }
  else {
    for(i=0; i<n_pop; i++) 
      cout << mf_rates[i] << " " ; 
    cout << endl ;
  } 
  
}

void open_files() {
  str_mean_rates = path + "/mean_rates.dat" ; 
  str_filter_rates = path + "/filter_rates.dat" ; 
  str_inputs = path + "/inputs.dat" ; 
  str_overlaps = path + "/overlaps.dat" ; 

  str_ff_inputs = path + "/ff_inputs.dat" ; 
  file_ff_inputs.open(str_ff_inputs.c_str(), ios::out | ios::ate) ; 
  
  file_mean_rates.open(str_mean_rates.c_str(), ios::out | ios::ate) ; 
  file_filter_rates.open(str_filter_rates.c_str(), ios::out | ios::ate) ; 
  file_inputs.open(str_inputs.c_str(), ios::out | ios::ate) ;   
  file_overlaps.open(str_overlaps.c_str(), ios::out | ios::ate) ; 
}

void close_files() { 
  file_mean_rates.close() ; 
  file_filter_rates.close() ; 
  file_inputs.close() ; 
  file_overlaps.close() ;
  file_ff_inputs.close() ;
}


void init_globals() { 
  m0 = M0 ; 
  M0_DT = m0 * DT ;
  
  n_neurons = (unsigned long) (n_neurons * 10000) ;     
  sqrt_K = sqrt( (float) K) ; 
  
  Ka = new float [n_pop]() ; 
  sqrt_Ka = new float [n_pop]() ;
  
  for(i=0;i<n_pop;i++) 
    if(IF_MATO_K) {
      Ka[i] = K * n_frac[i] ;
      sqrt_Ka[i] = sqrt( Ka[i] ) ;      
    }
    else {
      Ka[i] = K ; 
      sqrt_Ka[i] = sqrt_K ; 
    }
  
  n_per_pop = new unsigned long [n_pop]() ; 
  
  for(i=0;i<n_pop;i++) 
    n_per_pop[i] = (unsigned long) ( n_frac[i] * (float) n_neurons ) ; 
  
  cum_n_per_pop = new unsigned long [n_pop+1]() ; 
  
  cout <<"cum_n_per_pop=" << " " ; 
  for(i=0;i<n_pop+1;i++) { 
    for(j=0;j<i;j++) 
      cum_n_per_pop[i] += n_per_pop[j] ; 
    cout << cum_n_per_pop[i] << " " ;
  }
  cout << endl ; 
  
  Vr[0] = VrE ;
  if(n_pop==2)
    Vr[1] = VrI ;

  Vl[0] = VlE ;
  if(n_pop==2)
    Vl[1] = VlI ;
  
  A_CUE[0] = A_CUE_E ;
  if(n_pop==2)
    A_CUE[1] = A_CUE_I ;

  EPS_CUE[0] = EPS_CUE_E ;
  if(n_pop==2)
    EPS_CUE[1] = EPS_CUE_I ;

  A_ERASE[0] = A_ERASE_E ;
  if(n_pop==2)
    A_ERASE[1] = A_ERASE_I ;

  EPS_ERASE[0] = EPS_ERASE_E ;
  if(n_pop==2)
    EPS_ERASE[1] = EPS_ERASE_I ;
  
  
  which_pop = (int *) malloc( (unsigned long) n_neurons * sizeof(int) ) ; 
  
  for(i=0;i<n_pop;i++) 
    for(j=0; j<n_neurons; j++)
      if(j>=cum_n_per_pop[i] && j<cum_n_per_pop[i+1])
	which_pop[j] = i ; 
  
  mean_rates = new int [n_pop]() ; // population averaged rate also averaged over TIME_WINDOW  
  filter_rates = new int [n_neurons]() ; // temporal averaged over TIME_WINDOW

  if(REC_FF) {
    filter_ff_inputs = new float [n_neurons]() ; // temporal averaged over TIME_WINDOW
    for(j=0; j<n_neurons; j++)
      filter_ff_inputs[j] = 0 ;
  }
  // ISI = new float [n_neurons]() ; // temporal averaged over TIME_WINDOW 
  
  // h^(ab)_i=h^b_i, inputs from presynaptic population b to postsynaptic neuron (i,a)
  // here i goes threw all the populations 
  inputs = new float *[n_pop]() ; 
  for(i=0;i<n_pop;i++) // presynaptic population b 
    inputs[i] = new float [n_neurons]() ; 
  
  if(IF_NMDA) {
    inputs_nmda = new float *[n_pop]() ; 
    for(i=0;i<n_pop;i++) // presynaptic population b 
      inputs_nmda[i] = new float [n_neurons]() ; 
  }

  if(REC_INPUTS) {
    filter_inputs = new float *[n_pop]() ; 
    for(i=0;i<n_pop;i++) // presynaptic population b 
      filter_inputs[i] = new float [n_neurons]() ; 
  }
  
  // htot_i = h^E_i + h^I_i, net input into neuron i
  net_inputs = new float [n_neurons]() ; 
  if(IF_RK2)
    net_inputs_RK2 = new float [n_neurons]() ; 

  ff_inputs = new float [n_neurons]() ; 
  bg_inputs = new float [n_neurons]() ; 
  ext_inputs_scaled = new float [n_pop]() ; 
    
  mf_rates = new float [n_pop]() ; 
  mean_field_rates() ; 

  J_scaled = new float [n_pop*n_pop]() ;
  
  if(IF_NMDA)
    J_nmda = new float [n_pop*n_pop]() ;  
    
  if(IF_FF_LAYER) {
    J_FF = new float [n_neurons]() ; 
    J_FF0 = new float [N_POISSON]() ;
    J_task = new float [n_neurons]() ;
    
    if(IF_SPARSE_FF==0) 
      J_FF_all = new float [n_neurons*N_POISSON]() ;
  }

  if(IF_LOW_RANK) {
    overlaps = new float [n_pop]() ; 
    ksi = new float [n_per_pop[0]]() ;    
    if(RANK==2) {
      ksi_1 = new float [n_per_pop[0]]() ;
      overlaps_1 = new float [n_pop]() ;       
    }
    
    sample_A = new float [n_per_pop[0]]() ; 
    sample_B = new float [n_per_pop[0]]() ; 
    
    if(IF_SPEC) 
      if(RANK==2) {
	idx_perm_E = new unsigned long [n_per_pop[0]]() ; 
	idx_perm_I = new unsigned long [n_per_pop[1]]() ;
	idx_perm = new unsigned long [n_neurons]() ; 
	theta_1 = new float [n_neurons]() ; 
      }    
  }
  
  if(IF_GAUSS)
    prefactor = new float [n_pop*n_neurons]() ;   
  
}

void delete_globals() { 
  delete [] mean_rates ; 
  delete [] filter_rates ; 
  
  delete [] inputs ;
  if(IF_NMDA)
    delete [] inputs_nmda ;

  if(REC_INPUTS) delete [] filter_inputs ;
  
  delete [] net_inputs ;
  delete [] bg_inputs ;
  
  if(IF_RK2)
    delete [] net_inputs_RK2 ; 
  
  delete [] ext_inputs ;
  delete [] sigma_FF ; 
  delete [] ff_inputs ;
  
  delete [] J ; 
  delete [] J_scaled ;

  if(IF_NMDA)
    delete [] J_nmda ; 

  if(REC_FF)
    delete [] filter_ff_inputs ;
  
  if(IF_FF_LAYER) {    
    delete [] J_FF ;
    delete [] J_FF0 ;
    delete [] J_task ;
    
    if(IF_SPARSE_FF==0) 
      delete [] J_FF_all ;    
  }

  
  delete [] n_per_pop ;
  delete [] cum_n_per_pop ;
  free(which_pop) ;
  
  delete [] sqrt_Ka ;
  delete [] Ka ;
  
  if(IF_LOW_RANK) {
    delete [] ksi ; 
    delete [] ksi_init ;
    // delete [] ksi_scaled ;
    
    delete [] sample_A ; 
    delete [] sample_B ; 
    
    if(RANK==2) {
      delete [] ksi_1 ;
      // delete [] ksi_1_scaled ;
      delete [] overlaps_1 ;
      
    }
    delete [] overlaps ;
  }

  if(IF_GAUSS) 
    delete [] prefactor ; 

  if(IF_SPEC) {
    delete [] theta ;

    if(IF_FF_LAYER)
      delete [] theta_ff ;
    
    if(RANK==2) {
      // delete [] idx_perm ;
      // delete [] idx_perm_E ;
      // delete [] idx_perm_I ;
      delete [] theta_1 ;       
    }
  }
}

void make_dir(string path) {
  string mkdirp = "mkdir -p " ; 
  mkdirp += path ; 
  
  const char * cmd = mkdirp.c_str() ; 
  const int dir_err = system(cmd) ; 
  
  cout << path << endl ; 
}

#define VarToString(name) var_to_string(#name, (name))

template <class T>
const char* var_to_string(const char *name, const T& val) {
  // cout << name << " = " << val << endl;
  return name ;
}

void get_args(int argc , char** argv) { 
  if(argv[1] != NULL) { 
    n_pop = (int) atoi(argv[1]) ; 
    n_neurons = (unsigned long) atoi(argv[2]) ; 
    K = (float) atof(argv[3]) ;
    if(n_pop!=1)
      dir = argv[4] ;
  }
  else {
    cout << "n_pop ? " ;
    cin >> n_pop ;
    cout << "n_neurons ? " ;
    cin >> n_neurons ;
    cout << "K ? " ;
    cin >> K ;
    cout << "Directory ? " ;
    cin >> dir ;
  } 
}

///////////////////////////////////////////////////////////////////////

void get_param() {

  cout << "reading parameters from : " ;
  string file_name = "../parameters/" + to_string(n_pop) + "pop/" + dir +".txt" ; 
  cout << file_name << endl; 
  
  ext_inputs = new float [n_pop]() ;
  sigma_FF = new float [n_pop]() ;
  J = new float [n_pop*n_pop]() ;
  
  float *Tsyn ;
  Tsyn = new float [n_pop*n_pop]() ; 
  
  string token ;
  string::size_type sz;
  ifstream file(file_name.c_str()) ; 
  int i,j;
  
  i=0 ; 
  while(getline(file, token)) {
    j=-1 ; 
    istringstream line(token);
    while(line >> token) {
      if(i==0)
	if(j!=-1) ext_inputs[j] = stod(token, &sz) ; 
      
      if(i==1)
	if(j!=-1) J[j] = stod(token, &sz) ; 
      
      if(i==2)
	if(j!=-1) Tsyn[j] = stod(token, &sz) ; 
      
      if(i==3)
	if(j!=-1) sigma_FF[j] = stod(token, &sz) ; 
      
      j++ ; 
    }
    
    if(file.unget().get() == '\n') {
      i++ ;
    }
  }
  
  if(n_pop==1) { 
    J[0] = J0 ; 
    ext_inputs[0] = I0 ; 
  }
  
  cout << "ext_inputs" << endl ;
  for(int i=0;i<n_pop;i++)
    cout << ext_inputs[i] << " " ;
  cout << endl ;

  cout << "ext_var" << endl ;
  for(int i=0;i<n_pop;i++)
    cout << sigma_FF[i] << " " ;
  cout << endl ;
  
  cout << "J" << endl ;
  for(int i=0;i<n_pop;i++) {
    for(int j=0;j<n_pop;j++)
      cout << J[j+i*n_pop] << " " ;
    cout << endl ;
  }
  
  cout << "Tsyn" << endl ;
  for(int i=0;i<n_pop;i++) {
    for(int j=0;j<n_pop;j++)
      cout << Tsyn[j+i*n_pop] << " " ;
    cout << endl ;    
  }
  
  
}

float Phi(float x) { // Gaussian CDF
  return 0.5 * ( 1.0 + erf(x/sqrt(2.0) ) ) ; 
}

float threshold_linear(float x) {
  if(x>0.0) 
    return x ; 
  else 
    return 0. ; 
}

float cut_LR(float x) {
  if(x>1.0)
    x = 1.0 ;
  if(x<0.0)
    x = 0 ;
  
  return x ; 
}

void create_dir() { 
  
  path += "simulations/" ; 
  
  if(IF_LIF) 
    path += "lif/" ; 
  if(IF_BIN) 
    path += "binary/" ; 
    
  ostringstream str_I0, str_J0 ;
  str_I0 << fixed << setprecision(2) << I0 ;
  str_J0 << fixed << setprecision(2) << abs(J0) ;
  
  if(n_pop==1)
    dir = "I0_" + str_I0.str()  + "_J0_" + str_J0.str() ;
  
  path += to_string(n_pop)+"pop/" + dir  ;
  /* path += "/N" + to_string(n_neurons) ; */
  
  if(n_pop==1)
    path += "/N" + to_string(n_per_pop[0]/1000) ; 
  else 
    path += "/NE_" + to_string(n_per_pop[0]/1000) +  "_NI_" + to_string(n_per_pop[1]/1000) ; 
  
  path += "/K" + to_string((int)K) ; 

  ostringstream str_seed;  
  str_seed << SEED_CON ; 
  
  if(IF_CON_DIR)
    path += "/seed_" + str_seed.str() ;
  
  ostringstream str_tau_fac, str_tau_rec, str_use ;
  str_tau_fac << fixed << setprecision(0) << TAU_FAC[0] ; 
  str_tau_rec << fixed << setprecision(0) << TAU_REC[0] ; 
  str_use << fixed << setprecision(2) << USE[0] ; 

  if(IF_STP) 
    path += "/STP/Tf_" + str_tau_fac.str() + "_Tr_" +  str_tau_rec.str() + "_U_" +  str_use.str() ;  
  
  ostringstream str_kappa, str_kappa_1 ; 
  str_kappa << fixed << setprecision(2) << (float) KAPPA ; 
  str_kappa_1 << fixed << setprecision(2) << (float) KAPPA_1 ; 
  
  ostringstream str_map_seed ; 
  str_map_seed << fixed << setprecision(0) << MAP_SEED ; 

  ostringstream str_EE, str_EI, str_IE, str_II ;  
  str_EE << fixed << setprecision(0) << SIGMA[0] ; 
  str_EI << fixed << setprecision(0) << SIGMA[1] ; 
  str_IE << fixed << setprecision(0) << SIGMA[2] ; 
  str_II << fixed << setprecision(0) << SIGMA[3] ; 
  
  if(IF_STRUCTURE) { 
    // if(IF_SPEC) {
    //   if(RANK==1)
    // 	path += "/spec/kappa_" + str_kappa.str() ;
    //   if(RANK==2) {
    // 	path += "/spec/kappa_" + str_kappa.str() + "_kappa_1_" + str_kappa_1.str() ;
    // 	if(FIX_MAP_SEED) 
    // 	  path += "/seed_" + str_map_seed.str() ; 
    //   }
    // }
    
    if(IF_LOW_RANK) {
      path += "/low_rank/kappa_" + str_kappa.str() ; 
      if(RANK==2) 
	path += "_kappa_1_" + str_kappa_1.str() ;       
    }
    
    if(IF_RING)
      path += "/ring/kappa_" + str_kappa.str() ;
    
    if(IF_GAUSS)
      path += "/gauss/EE_" + str_EE.str() + "_EI_" + str_EI.str() +"_IE_" + str_IE.str() + "_II_" + str_II.str() ;     
  } 
  
  ostringstream str_kappa_ext, str_phi_ext, str_phi_dist ; 
  str_kappa_ext << fixed << setprecision(2) << KAPPA_EXT ; 
  str_phi_ext << fixed << setprecision(3) << PHI_EXT ; 
  
  // christos 
  ostringstream str_phi_cue, str_phi_erase ; 
  str_phi_cue << fixed << setprecision(3) << PHI_CUE ; 
  str_phi_erase << fixed << setprecision(3) << PHI_ERASE ; 
  
  ostringstream str_A_cue, str_eps_cue ; 
  str_A_cue << fixed << setprecision(2) << A_CUE[0] ; 
  str_eps_cue << fixed << setprecision(2) << EPS_CUE[0] ; 

  ostringstream str_A_erase, str_eps_erase ; 
  str_A_erase << fixed << setprecision(2) << A_ERASE[0] ; 
  str_eps_erase << fixed << setprecision(2) << EPS_ERASE[0] ; 

  if(IF_DPA || IF_DUAL || IF_DRT) {
    if(IF_DPA) 
      path += "/DPA/kappa_" + str_kappa_ext.str() ; 
    if(IF_DUAL) 
      path += "/dual_task/kappa_" + str_kappa_ext.str() ; 
    if(IF_DRT) 
      path += "/DRT/kappa_" + str_kappa_ext.str() ;

    if(SAMPLE)
      path += "/sample_B" ;
    else
      path += "/sample_A" ;
    
    if(DISTRACTOR && (IF_DUAL || IF_DRT) )
      path += "/NoGo" ; 
    else 
      path += "/Go" ; 
  }
  
  if(IF_CHRISTOS){
    path += "/christos/cue_A_" +  str_A_cue.str() + "_eps_" + str_eps_cue.str() + "_phi_" + str_phi_cue.str() ;
    path += "/dist_A_" +  str_A_erase.str() + "_eps_" + str_eps_erase.str() + "_phi_" + str_phi_erase.str() ; 
  }
  

  ostringstream str_m0 ; 
  str_m0 << fixed << setprecision(4) << M0 ; 
  if(IF_LOOP_M0) 
    path += "/m0_" + str_m0.str() ; 

  ostringstream str_gain ; 
  str_gain << fixed << setprecision(1) << GAIN ; 
  if(IF_LOOP_GAIN) 
    path += "/gain_" + str_gain.str() ; 
  
  if(IF_TRIALS) 
    path += "/trial_" + to_string( (int) TRIAL_ID ) ;
  
  if(IF_INI_COND) 
    path += "/ini_cond_" + to_string( (int) INI_COND_ID ) ; 
  
  make_dir(path) ; 
  
  cout << "Created directory : " ; 
  cout << path << endl ; 

}

template <class T> 
void read_from_file(string path, string file_name, T * &array, size_t array_size) {

  int dum ; 
  
  string file_path ; 
  file_path = path + "/" + file_name + ".dat" ; 
  cout << "reading from: " << file_path << endl ; 
  
  struct stat buffer ; 
  FILE *file ; 
  
  if (stat (file_path.c_str(), &buffer) == 0) { 
    file = fopen(file_path.c_str(), "rb") ; 
    dum = fread(&array[0], sizeof array[0], array_size, file) ; 
    fclose(file) ; 
  } 
  else { 
    cout << "ERROR: " << file_name << ".dat not found" << endl ; 
    exit(-1) ; 
  }
}

template <class T> 
void write_to_file(string path, string file_name, T * &array, size_t array_size) {

  int dum ; 
  
  string file_path ; 
  file_path = path + "/" + file_name + ".dat" ; 
  cout << "writing to: " << file_path << endl ; 
  
  struct stat buffer ; 
  FILE *file ; 
  
  file = fopen(file_path.c_str(), "wb") ; 
  dum = fwrite(&array[0], sizeof array[0], array_size, file) ; 
  fclose(file) ;
  
}

void save_to_file() { 

  file_mean_rates << t_time - TIME_STEADY ; 
    
  for(i=0;i<n_pop;i++) { 
    file_mean_rates << " " << ( (float) mean_rates[i] ) * 1000. / TIME_WINDOW / (float) n_per_pop[i] ; 
    mean_rates[i] = 0 ; 
  }
  
  file_mean_rates << endl ; 
  
  if(IF_LOW_RANK) {
    file_overlaps << t_time - TIME_STEADY ; 
    for(i=0;i<n_pop;i++) { 
      file_overlaps << " " << overlaps[i] * 1000. / TIME_WINDOW * IS_STRUCT_SYN[i] ; 
      /* file_overlaps << " " << overlaps[i]*1000./TIME_WINDOW/(float)n_per_pop[i] ;  */
      overlaps[i] = 0 ; 
    } 
    file_overlaps << endl ; 
  }
  
  // filtered rates over tw
  file_filter_rates << t_time - TIME_STEADY ; 
  
  for(i=0;i<n_neurons;i++) { 
    file_filter_rates << " " << ( (float) filter_rates[i] ) * 1000. / TIME_WINDOW ; 
    filter_rates[i] = 0 ; 
  } 
  file_filter_rates << endl ; 
  
  // filtered inputs over tw 
  file_inputs << t_time - TIME_STEADY ;
  
  if(REC_INPUTS){
    for(i=0;i<n_pop;i++) 
      for(j=0;j<n_neurons;j++) { 
	file_inputs << " " << filter_inputs[i][j]*DT/TIME_WINDOW ; 
	filter_inputs[i][j] = 0 ; 
      } 
    file_inputs << endl ;
  }
}

void get_m1_phase() { 
  
  float dPhi = 0 ; 
  float xCord = 0, yCord = 0 ; 
  
  for(int i_pop=0; i_pop<n_pop; i_pop++) {
    
    dPhi = M_PI / (float) n_per_pop[i_pop] ; 
    xCord = 0;
    yCord = 0 ; 
    
    for(unsigned long j=cum_n_per_pop[i_pop]; j < cum_n_per_pop[i_pop+1]; j++) {
      xCord += (float) filter_rates[j] * cos(2.0 * j * dPhi) / TIME_WINDOW * 1000.0 ; 
      yCord += (float) filter_rates[j] * sin(2.0 * j * dPhi) / TIME_WINDOW * 1000.0 ; 
    }
    
    m1[i_pop] = ( 2.0 / (float) n_per_pop[i_pop]) * sqrt(xCord * xCord + yCord * yCord) ; 
    phase[i_pop] = 0.5 * atan2(yCord, xCord) ; 
    
    if(phase[i_pop] < 0)
      phase[i_pop] = phase[i_pop] + M_PI ;
    
    // phase[i_pop] *= 180.0/M_PI ; 
    phase[i_pop] *= 180.0/M_PI - 180. ; 
  }
}

#endif
