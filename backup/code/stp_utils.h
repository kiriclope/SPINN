#ifndef __STPUTILS__ 
#define __STPUTILS__ 

//STP Parameters 

float *u_stp ; // availability variable 
float *x_stp, *x_stp_old ; // resource variable 
float *A_u_x_stp ; // A*u*x variable 

float *u_stp_ff ; // availability variable 
float *x_stp_ff ; // resource variable 
float *A_u_x_stp_ff ; // A*u*x variable 

float *u_stp_mongillo, *x_stp_mongillo ; // availability variable 

float *unbind_rate, *refill_rate, *bind_rate ;

string str_u_stp ; 
ofstream file_u_stp ; 

string str_x_stp ; 
ofstream file_x_stp ; 

void init_stp_globals(){ 
  u_stp = new float [n_neurons]() ; 
  x_stp = new float [n_neurons]() ; 
  A_u_x_stp = new float [n_neurons]() ; 
  
  if(IF_STP_FF) {
    u_stp_ff = new float [N_POISSON]() ; 
    x_stp_ff = new float [N_POISSON]() ; 
    A_u_x_stp_ff = new float [n_neurons]() ; 
  }
  
  /* u_stp_mongillo = new float [n_neurons*n_neurons]() ;  */
  /* x_stp_mongillo = new float [n_neurons*n_neurons]() ;  */
  
  /* x_stp_old = new float [n_neurons]() ;  */
  
  /* unbind_rate = new float [n_neurons]() ;  */
  /* refill_rate = new float [n_neurons]() ;  */
  /* bind_rate = new float [n_neurons]() ;  */
  
}

void delete_stp_globals(){ 
  delete [] u_stp ; 
  delete [] x_stp ;

  if(IF_STP_FF){
    delete [] u_stp_ff ; 
    delete [] x_stp_ff ;
    delete [] A_u_x_stp_ff ;
  }
  
  /* delete [] x_stp_old ;  */
  
  /* delete [] u_stp_mongillo ;  */
  /* delete [] x_stp_mongillo ; */
  
  /* delete [] unbind_rate ; */
  /* delete [] refill_rate ; */
  /* delete [] bind_rate ;  */
}

void open_stp_files(){ 
  str_u_stp = path + "/u_stp.dat" ; 
  file_u_stp.open(str_u_stp.c_str(), ios::out | ios::ate) ;
  
  str_x_stp = path + "/x_stp.dat" ; 
  file_x_stp.open(str_x_stp.c_str(), ios::out | ios::ate) ; 
}

void close_stp_files(){
  file_u_stp.close() ; 
  file_x_stp.close() ; 
}

void markram() {
  
  x_stp[i_neuron] -= u_stp[i_neuron] * x_stp[i_neuron] ; 
  u_stp[i_neuron] += USE[pre_pop] * (1.0 - u_stp[i_neuron]) ; 
  A_u_x_stp[i_neuron] = u_stp[i_neuron] * x_stp[i_neuron] ;
  
}

/* Markram98, Mato  */
void mato() { 

  if(TAU_FAC[pre_pop]!=0) 
    u_stp[i_neuron] = u_stp[i_neuron] * exp(- ISI / TAU_FAC[pre_pop] ) + USE[pre_pop] * ( 1. - u_stp[i_neuron] * exp(- ISI / TAU_FAC[pre_pop] ) ) ; 
  
  x_stp[i_neuron] = x_stp[i_neuron] * ( 1. - u_stp[i_neuron] ) * exp(- ISI / TAU_REC[pre_pop] ) + 1. - exp(- ISI / TAU_REC[pre_pop] ) ; 
  A_u_x_stp[i_neuron] = u_stp[i_neuron] * x_stp[i_neuron] ; 
}

void mato_ff() {
  if(TAU_FAC_FF!=0)
    u_stp_ff[i_neuron_ff] = u_stp_ff[i_neuron_ff] * exp(- ISI_FF / TAU_FAC_FF ) + USE_FF * ( 1. - u_stp_ff[i_neuron_ff] * exp(- ISI_FF / TAU_FAC_FF ) ) ;
  
  x_stp_ff[i_neuron_ff] = x_stp_ff[i_neuron_ff] * ( 1. - u_stp_ff[i_neuron_ff] ) * exp(- ISI_FF / TAU_REC_FF ) + 1. - exp(- ISI_FF / TAU_REC_FF ) ; 
  A_u_x_stp_ff[i_neuron_ff] = u_stp_ff[i_neuron_ff] * x_stp_ff[i_neuron_ff] ; 
}

/* void release_stp() { */
  
/*   if(IF_MARKRAM) {  */
/*     u_stp[i_neuron] *= EXP_DT_TAU_FAC[pre_pop] ;  */
/*     /\* x_stp_old[i_neuron] = x_stp[i_neuron] ;  *\/ */
/*     x_stp[i_neuron] += DT_OVER_TAU_REC[pre_pop] * (1.0 - x_stp[i_neuron]) ;  */
/*   } */
  
/*   if(IF_MONGILLO==1) { */
/*     if( u_stp[i_neuron] == 1.0 && unif(rand_gen) < DT / TAU_FAC[pre_pop]) {// u = 1 -> 0 , calcium unbinds with rate 1/Tfac  */
/*       u_stp[i_neuron] = 0.0 ;  */
/*       /\* unbind_rate[i_neuron] += 1.0 ;  *\/ */
/*     } */
/*     if( x_stp[i_neuron] == 0.0 && unif(rand_gen) < DT / TAU_REC[pre_pop]) {// x = 0 -> 1 , neurotransmitter refills with rate 1/Trec  */
/*       x_stp[i_neuron] = 1.0 ; */
/*       /\* refill_rate[i_neuron] += 1.0 ;  *\/ */
/*     } */
/*   } */
/* } */

/* void mongillo() { // mongillo  */
  
/*   // Updating STP variables */
/*   if( u_stp[i_neuron]==0.0 && unif(rand_gen) < USE[pre_pop]*DT) {  */
/*     u_stp[i_neuron] = 1.0 ; // calcium binds with probability Use u = 0->1 */
/*     /\* bind_rate[i_neuron] += 1.0 ;  *\/ */
    
/*     if(x_stp[i_neuron]==1.0) { */
/*       A_u_x_stp[i_neuron] = 1.0 ; x_stp[i_neuron] = 0.0 ; */
/*     } // neurotransmitter release if available x = 1->0  */
/*     else */
/*       A_u_x_stp[i_neuron] = 0.0 ; */
/*   } */
/*   else */
/*     A_u_x_stp[i_neuron] = 0.0 ; */
/* } */

/* void mongillo_alt() { // mongillo */
  
/*   u_stp[i_neuron] = u_stp_mongillo[i_neuron+id_post[j]*n_neurons] ;  */
/*   x_stp[i_neuron] = x_stp_mongillo[i_neuron+id_post[j]*n_neurons] ;  */
  
/*   if (u_stp[i_neuron]==1.0) */
/*     if ( unif(rand_gen)<(1.0-exp(-ISI/TAU_FAC[pre_pop]))*DT) */
/*       u_stp[i_neuron]=0.0; */
  
/*   if (u_stp[i_neuron]==0.0) */
/*     if (unif(rand_gen)<USE[pre_pop]*DT) */
/*       u_stp[i_neuron]=1.0; */
  
/*   if (x_stp[i_neuron]==0.0) */
/*     if ( unif(rand_gen)<(1.0-exp(-ISI/TAU_REC[pre_pop]))*DT) */
/*       x_stp[i_neuron]=1.0 ;  */
  
/*   A_u_x_stp[i_neuron] = u_stp[i_neuron] * x_stp[i_neuron] ; */
  
/*   if (u_stp[i_neuron]==1.0) */
/*     u_stp_mongillo[i_neuron+id_post[j]*n_neurons] = 1.0 ; */
/*   else */
/*     x_stp_mongillo[i_neuron+id_post[j]*n_neurons] = x_stp[i_neuron] ;  */
/* } */

void update_stp_variables_lif() {
  
  if(IF_MARKRAM)
    markram() ; 
  if(IF_MATO)
    mato() ; 
  /* if(IF_MONGILLO==1) */
  /*   mongillo() ; */
  /* if(IF_MONGILLO==2) */
  /*   mongillo_alt() ;  */
  
}

void save_xy_to_file() { 
  
  /* cout << "unbind_rate" ;  */
  /* for(j=0;j<5;j++) {  */
  /*   cout << unbind_rate[j] *1000./TIME_WINDOW << " " ;  */
  /*   unbind_rate[j] = 0 ;  */
  /* } */
  /* cout << endl ;  */
  
  /* cout << "refill rate" ; */
  /* for(j=0;j<5;j++) { */
  /*   cout << refill_rate[j] *1000./TIME_WINDOW << " " ; */
  /*   refill_rate[j] = 0 ;     */
  /* } */
  /* cout << endl ;  */
  
  /* cout << "bind rate" ;  */
  /* for(j=0;j<5;j++) {  */
  /*   cout << bind_rate[j] *1000./TIME_WINDOW << " " ;  */
  /*   refill_rate[j] = 0 ;  */
  /* } */
  /* cout << endl ;  */

  file_u_stp << t_time-TIME_STEADY ; 
  for(j=0;j<n_neurons;j++) 
    file_u_stp << " " << u_stp[j] ;  
  file_u_stp << endl ; 
  
  file_x_stp << t_time-TIME_STEADY ; 
  for(j=0;j<n_neurons;j++) 
    file_x_stp << " " << x_stp[j] ;  
  file_x_stp << endl ; 
}


#endif
