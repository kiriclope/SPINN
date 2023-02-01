#ifndef __GLOBALVARS__
#define __GLOBALVARS__

using namespace:: std ;
using namespace libconfig;

random_device rd ;
mt19937 rand_gen(rd()) ; // random device

uniform_real_distribution<float> unif(0.0, 1.0) ;
normal_distribution<float> white_noise(0.0, 1.0) ;

#define IF_BIN 0
// binary network
#define IF_LIF 1 
// LIF network
#define IF_RATE 0
// Rate model

#define E_frac (float) 0.8 
// fraction in excitatory pop
const float n_frac[2] = { E_frac, round( (1.0 - E_frac)*100.0) / 100.0 } ; 
// fraction of neurons in each pop

#define IF_MATO_K 1 
// scaling as in Mato et al. 
#define IF_RESCALE 0 
// rescaling by multiplying Jijs by Vth-Vl 

////////////////////////////////////////////////////////////////////////////////////
// Simulation globals 
////////////////////////////////////////////////////////////////////////////////////
#define DT (float) .1 

#define TIME_INI (double) 0E3 
#define TIME_STEADY (double) 5E3 

#define DURATION (double) 30E3 

#define TIME_WINDOW (double) 2.0E3 
#define TIME_REC (double) 100E3 
#define TIME_REC_SPIKES (double) 0E3 

#define REC_INPUTS 0
#define IF_RK2 0 

#define N_PREF 10000 

#define IF_TRIALS 1 
int TRIAL_ID = 0 ;

#define IF_INI_COND 0 
/* #define INI_COND_ID 0  */
int INI_COND_ID = 0 ;

#define IF_SAVE_VOLT 0 

#define IF_CON_DIR 0 
#define SEED_CON 1 
mt19937 con_gen(sqrt(SEED_CON)) ; // random device for gen connectivity

////////////////////////////////////////////////////////////////////////////////////
// Network globals 
////////////////////////////////////////////////////////////////////////////////////
const float eps[2] = {1.0,-1.0} ; 
const float Trate[2] = {10.0, 10.0} ; 

#define J0 -1.0 // for one population 
#define I0 1.0 // for one population 
#define Tsyn0 1.0 // for one population

#define GAIN (float) 1.0 
float M0 = 0.004 ; 
float M0_DT=0 ; 

#define IF_LOOP_M0 0 
#define IF_LOOP_GAIN 0 

float m0, m1[2], phase[2] ; 
float duration, time_rec, time_rec_spikes, time_steady ; 

string dir ; 
string path = "../" ; 

unsigned long i, j, k, l, i_neuron, i_neuron_ff ; 
float t_time, percentage, t_window=0. ; 
int pre_pop, post_pop ; 

int n_pop ; 
unsigned long n_neurons ; 
float K, sqrt_K, *sqrt_Ka, *Ka ; 

float *ext_inputs, *ff_inputs, *bg_inputs, *ext_inputs_scaled ; 
float *J, *J_scaled, *J_nmda ; 

unsigned long *n_per_pop, *cum_n_per_pop ; 
int *which_pop ; 

float *mf_rates ; 

float *volt, vold ;
float RK1, RK2 ; 
float *t_spike ; 
float *t_spike_ff ; 

int *mean_rates, *filter_rates ; 

float **inputs, **inputs_nmda ; 
float **filter_inputs ; 
float *filter_ff_inputs ;

float *net_inputs ; 
float *net_inputs_RK2 ; 

///////////////////
// External inputs 
///////////////////
int IF_NOISE = 0 ;

#define IF_GAUSS_NOISE 0
#define IF_SQRT_K_NOISE 0
#define IF_COS_NOISE 0
double noisy_phase ; 

#define IF_NO_MAP 0 
#define IF_TUNED_EXT 0


////////////////////////////////////////////////////////////////////////////////////
// Recurrent connectivity
////////////////////////////////////////////////////////////////////////////////////

#define IF_GEN_CON 0 
#define IF_SAVE_CON_VEC 0 
#define IF_SAVE_SPARSE_REP 0 
#define IF_CHECK_SPARSE_REP 0 

string con_path = "../connectivity/" ;

int *n_post, *avg_n_post, *con_vec;
int **n_pre, *avg_n_pre ; 
unsigned long *id_post, *idx_post ; 
unsigned long total_n_post = 0 ; 

float *K_over_Na ; 
float *con_prob ; 
float *prefactor ; 

int IF_STRUCTURE ; 
const float IS_STRUCT_SYN[4] = {1.0, 0.0, 0.0, 0.0} ; 
float *theta, *theta_1, *theta_ff ; 

///////////
//Ring
///////////
#define IF_RING 0
// same as van Wreeswijk, Les Houches 
#define IF_SPEC 0 
// same as Rao et al.
float KAPPA = 4.5 ;
// 4.5 for chengyu's task
// 0.25 for christos
// first fourier of the connection proba
// 4 // 3.5 // 12.0 14 // rank 1: 4.0 
const float kappas[4] = {KAPPA, (float) KAPPA * (float) 0.25 ,  KAPPA,  (float) KAPPA * (float) 0.25 } ;
float kappa, kappa_K_N ;  

///////////
//Gauss
///////////
#define IF_GAUSS 0
float *X ; 
const float SIGMA[4] = {60.0, 60.0, 70.0, 60.0} ; 
#define DEG_TO_RAD (float) M_PI/180.0 
/* #define L 2.0*M_PI  */

/////////////
// low-rank
/////////////
#define IF_LOW_RANK 1 
#define RANK 2 
#define IF_GEN_KSI 0 

float *overlaps, *overlaps_1 ; 

#define FIX_KSI_SEED 1 
#define SEED_KSI (float) 2 
string ksi_path ; 

float KAPPA_1 = 4.5 ;
#define MAP_ANGLE (float) M_PI/4.0 
// in radians 

#define FIX_MAP_SEED 1 
#define MAP_SEED 4 

mt19937 ksi_gen(exp(SEED_KSI)) ; 
mt19937 ksi_1_gen(sqrt(SEED_KSI)) ; 
mt19937 covar_ksi_gen(10.0) ; 

// // cov for last days
// float COV_MAT[16] = { // sample A, sample B, ksi, ksi1
// 		     1.0, 0.0, 0.5, 0.0,
// 		     0.0, 1.0, -0.5, 0.0,
// 		     0.5, -0.5, 1.0, 0.0,
// 		     0.0, 0.0, 0.0, 1.0,
// } ;

/* // cov for early days  */
float COV_MAT[16] = { // sample A, sample B, ksi, ksi1
		     1.0, 0.0, 0.6, 0.2,
		     0.0, 1.0, -0.6, 0.2,
		     0.6, -0.6, 1.0, -0.1 / 2000.0,
		     0.2, 0.2, -0.1 / 2000.0, 1.0,
} ;

float *ksi, *ksi_1, *ksi_init, *ksi_scaled, *ksi_1_scaled ; 
float *sample_A, *sample_B ; 
unsigned long *idx_perm, *idx_perm_E, *idx_perm_I ; 
float kappa_1, kappa_1_K_N ; 

////////////////////////////////////////////////////////////////////////////////////
// feedforward layer
////////////////////////////////////////////////////////////////////////////////////
#define IF_FF_LAYER 0
unsigned long *id_post_ff, *idx_post_ff ; 
int *n_post_ff, *avg_n_post_ff, *con_vec_ff ; 
unsigned long total_n_post_ff = 0 ;
float *con_prob_ff ; 
float K_over_N_FF ; 
float *prefactor_ff ; 

#define IF_FF_EI 0
#define IF_SPARSE_FF 0
#define IF_SYN_DYN_FF 0 // 1 exponential synapses, 0 delta synapses

string con_path_ff = "../connectivity/" ; 

#define IF_GEN_CON_FF 0 
#define IF_SAVE_SPARSE_REP_FF 0 

#define N_POISSON 32000 
#define REC_FF 0

#define K_FF (float) 2000 
#define GAIN_FF (float) 0.25

float NU_FF = 0 ;

float *J_FF, *J_FF0, *J_FF_all, *J_task ; 
int poisson_rates=0 ; 

float *sigma_FF ; 
float var_ff[2] = {0.0, 0.0} ; 

///////////
//Ring
///////////
#define IF_RING_FF 0
#define IF_SPEC_FF 0 

#define KAPPA_FF (float) .125 

///////////
// Gauss
///////////
#define IF_GAUSS_FF 0
#define IF_GAUSS_SPEC_FF 0

////////////////////////////////////////////////////////////////////////////////////
// LIF globals 
////////////////////////////////////////////////////////////////////////////////////

#define IF_SYN_DYN 1 
// 1 exponential synapses, 0 delta synapses 

float Vl[2] ; 
float VlE = 0.0 ;
#define VlI (float) 0.0 

float Vr[2] ; 
float VrE = -3.33 ; 
#define VrI -3.33 

#define Vth 20 
#define Vpeak 20. 
// Spikes Peak

float ISI ; 

string str_volt ; 
ofstream file_volt ; 

string str_spike_times ; 
ofstream file_spike_times ; 

const float TAU_MEM[2] = {20.0, 10.0} ; 
const float dt_over_tau_mem[2] = { (float) DT/TAU_MEM[0], (float) DT/TAU_MEM[1] } ; 
const float one_minus_dt_over_tau_mem[2]={ (float) (1.0-dt_over_tau_mem[0]), (float) (1.0-dt_over_tau_mem[1] )} ; 
const float one_plus_dt_over_tau_mem[2]={ (float) (1.0+dt_over_tau_mem[0]), float(1.0+dt_over_tau_mem[1]) } ; 
const float dt_over_dt_tau_mem[2]={(float) DT/ (float) (1.0+dt_over_tau_mem[0]), (float) DT/ (float) (1.0+dt_over_tau_mem[1]) } ; 

const float EXP_DT_TAU_MEM[2] = { (float) exp(-DT/TAU_MEM[0]) , (float) exp(-DT/TAU_MEM[1]) } ; 

const float TAU_SYN[4] = {4.0, 2.0, 4.0, 2.0} ; 
const float EXP_DT_TAU_SYN[4] = { (float) exp(-DT/TAU_SYN[0]) , (float) exp(-DT/TAU_SYN[1]), (float) exp(-DT/TAU_SYN[2]), (float) exp(-DT/TAU_SYN[3])} ; 

////////////
// NMDA
////////////
#define IF_NMDA 1 
const float TAU_NMDA[4] = {80.0, 40.0, 80.0, 40.0} ; 
const float EXP_DT_TAU_NMDA[4] = { (float) exp(-DT/TAU_NMDA[0]) , (float) exp(-DT/TAU_NMDA[1]), (float) exp(-DT/TAU_NMDA[2]), (float) exp(-DT/TAU_NMDA[3])} ; 
const float R_NMDA[2] = {1.0, 9.0} ; 

////////////////////////////////////////////////////////////////////////////////////
// STP globals 
////////////////////////////////////////////////////////////////////////////////////
#define IF_STP 1 

#define IF_MARKRAM 0
#define IF_MATO 1 
#define IF_MONGILLO 0 

const int stp_synapse[4] = {1, 0, 0, 0} ; 

const float TAU_FAC[2] = {400.0, 0.0} ; 

const float EXP_DT_TAU_FAC[2] = { (float) exp(-DT/TAU_FAC[0]) ,
				  (float) exp(-DT/TAU_FAC[1]) } ; 

const float TAU_REC[2] = {200.0, 800.0} ; 
  
const float DT_OVER_TAU_REC[2] = { (float) DT/TAU_REC[0] ,
				   (float) DT/TAU_REC[1] } ; 

const float USE[2] = {0.03, 0.5} ; 


#define IF_STP_FF 0
const float TAU_FAC_FF= 0 ; 
const float TAU_REC_FF= 800 ; 
const float USE_FF= 0.5 ; 
float ISI_FF=0 ; 

////////////////////////////////////////////////////////////////////////////////////
// Stimulus / Task globals 
////////////////////////////////////////////////////////////////////////////////////

int IF_STIM = 1 ; 

int SWITCH_ON = 0 ; 
int SWITCH_OFF = 0 ;

#define IF_TRACK 0 

#define KAPPA_EXT (float) 3.0 
#define PHI_EXT (float) 0.25 

#define KAPPA_DIST (float) 0.75 // 0.25 Go
#define PHI_DIST (float) 0.75 

#define KAPPA_CUE (float) 0.5 // 0.25 Go
#define KAPPA_TEST (float) 0.25 

#define IF_STEP 1 
#define T_STEP_ON (float) 2000.0
#define T_STEP_OFF (float) 3000.0
const float A_STEP[2] = {3.0, 0.0} ; 

////////////////////////////////////////////////////////////////////////////////////
// Christos
////////////////////////////////////////////////////////////////////////////////////

#define CUE (float) 1.0

float IF_CHRISTOS = 0 ; 
#define T_CUE_ON (float) 2000 
#define T_CUE_OFF (float) 3000 

#define T_ERASE_ON (float) 5000
#define T_ERASE_OFF (float) 6000 

float PHI_CUE = 0.25 ;

#define A_CUE_E 0.25 
#define A_CUE_I 0.0 
float A_CUE[2] ;

#define EPS_CUE_E 0.25 
#define EPS_CUE_I 0.0 
float EPS_CUE[2] ;

#define IF_DIST 0 
float PHI_ERASE = 0.75 ; 

#define A_ERASE_E 0.0 
#define A_ERASE_I 0.0 
float A_ERASE[2] ;

#define EPS_ERASE_E 0.25 
#define EPS_ERASE_I 0.0 
float EPS_ERASE[2] ;

#define IF_SECOND 0 

////////////////////////////////////////////////////////////////////////////////////
// DUAL TASK
////////////////////////////////////////////////////////////////////////////////////

float IF_DPA = 0 ;
float IF_DUAL = 0 ;
float IF_DRT = 0 ;

int SAMPLE = 0 ;
int DISTRACTOR = 0 ;

float EPS[2] = {1, -1} ;

float kappa_ext, kappa_dist, kappa_cue, kappa_test ; 

#define T_SAMPLE_ON (float) 2000 
#define T_SAMPLE_OFF (float) 3000 

#define T_DIST_ON (float) 4500 //4500 
#define T_DIST_OFF (float) 5500 // 5500 

#define T_RWD_ON (float) 6500 
#define T_RWD_OFF (float) 7500 

#define T_TEST_ON (float) 9000 
#define T_TEST_OFF (float) 10000 

#endif 
