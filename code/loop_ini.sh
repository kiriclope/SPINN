#!/usr/bin/env bash 

rand=$RANDOM # has to be outside

./bash_utils.sh $rand

temp_main=$(printf 'temp_%d_main.cpp' $rand) 
temp_globals=$(printf 'temp_%d_globals.h' $rand) 
temp_out=$(printf 'temp_%d_a' $rand) 

echo $temp_globals $temp_main $temp_out

IF_LIF=1 
IF_BIN=0 

IF_STP=1 
IF_GEN_CON=0 

RANK=1 
IF_GAUSS=0
IF_RING=1

IF_SPEC=0
FIX_MAP_SEED=1 

IF_LOW_RANK=0 
FIX_KSI_SEED=1 

sed -ie "s/ DURATION .*/ DURATION (double) 10E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_INI .*/ TIME_INI (double) 0E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_STEADY .*/ TIME_STEADY (double) 6E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_WINDOW .*/ TIME_WINDOW (double) .250E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_REC .*/ TIME_REC (double) 60E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_REC_SPIKES .*/ TIME_REC_SPIKES (double) 10E3 /" "$temp_globals" ; 

sed -ie "s/ IF_LIF .*/ IF_LIF ${IF_LIF} /" "$temp_globals" ; 
sed -ie "s/ IF_BIN .*/ IF_BIN ${IF_BIN} /" "$temp_globals" ; 

sed -ie "s/ IF_STP .*/ IF_STP ${IF_STP} /" "$temp_globals" ; 

sed -ie "s/ IF_TRIALS .*/ IF_TRIALS 0 /" "$temp_globals" ; 
sed -ie "s/ IF_INI_COND .*/ IF_INI_COND 1 /" "$temp_globals" ; 

sed -ie "s/ IF_GEN_CON .*/ IF_GEN_CON ${IF_GEN_CON} /" "$temp_globals" ; 
sed -ie "s/ IF_SAVE_CON_VEC .*/ IF_SAVE_CON_VEC 0 /" "$temp_globals" ; 
sed -ie "s/ IF_SAVE_SPARSE_REP .*/ IF_SAVE_SPARSE_REP 0 /" "$temp_globals" ; 

sed -ie "s/ RANK .*/ RANK ${RANK} /" "$temp_globals" ; 

sed -ie "s/ IF_SPEC .*/ IF_SPEC ${IF_SPEC} /" "$temp_globals" ; 
sed -ie "s/ FIX_MAP_SEED .*/ FIX_MAP_SEED ${FIX_MAP_SEED} /" "$temp_globals" ; 

sed -ie "s/ IF_LOW_RANK .*/ IF_LOW_RANK ${IF_LOW_RANK} /" "$temp_globals" ; 
sed -ie "s/ FIX_KSI_SEED .*/ FIX_KSI_SEED ${FIX_KSI_SEED} /" "$temp_globals" ; 

sed -ie "s/ IF_HYSTERESIS .*/ IF_HYSTERESIS 0 /" "$temp_globals" ; 
sed -ie "s/ IF_DPA .*/ IF_DPA 0 /" "$temp_globals" ; 
sed -ie "s/ IF_DUAL .*/ IF_DUAL 0 /" "$temp_globals" ; 

sed -ie "s/ IF_CHRISTOS .*/ IF_CHRISTOS 1 /" "$temp_globals" ; 
sed -ie "s/ IF_STEP .*/ IF_STEP 0 /" "$temp_globals" ; 

sed -ie "s/ IF_CON_DIR .*/ IF_CON_DIR 0 /" "$temp_globals" ; 

read n_pop N K dir n_ini <<<  "$1 $2 $3 $4 $5" 

sed -ie "s/ A_CUE_E .*/ A_CUE_E .25 /" "$temp_globals"  
sed -ie "s/ EPS_CUE_E .*/ EPS_CUE_E .2 /" "$temp_globals"  
    
sed -ie "s/ PHI_CUE (double) .*/ PHI_CUE (double) .25 /" "$temp_globals" ; 
sed -ie "s/ PHI_ERASE (double) .*/ PHI_ERASE (double) .75 /" "$temp_globals" ; 

for ini in $(seq 1 1 $n_ini); do
    
    echo "#########################################################################" 
    ./mem_usage.sh 
    # ./cpu_usage.sh 
    echo "#########################################################################" 
       
    echo "#########################################################################" 
    echo "simulation parameters:" 
    echo "n_pop ${n_pop} n_neurons ${N}0000 K ${K} ${dir} ini ${ini}" 
    echo "#########################################################################" 
    
    sed -ie "s/ INI_COND_ID .*/ INI_COND_ID ${ini} /" "$temp_globals" 
    
    echo "#########################################################################" 
    echo "compiling far condition"
    echo "#########################################################################"
    
    g++ -L/home/leon/bebopalula/cpp/libs/gsl/lib -I/home/leon/bebopalula/cpp/libs/gsl/include -std=c++11 ${temp_main} -Ofast -s -o ${temp_out}.out -lgsl -lgslcblas
	
    screen -dmS off_${n_pop}_pop_${dir}_N_${N}_K_${K}_ini_${ini}_off_far ./${temp_out}.out $n_pop $N $K ${dir}_off 
    
    ./mem_usage.sh 
    screen -dmS off_noise_0.75${n_pop}_pop_${dir}_N_${N}_K_${K}_ini_${ini}_off_far ./${temp_out}.out $n_pop $N $K ${dir}_noise_0.75_off 
    
    # sed -ie "s/ SIGMA_FF .*/ SIGMA_FF (double) 0.5 /" "$temp_globals" ;     
    # g++ -L/home/leon/bebopalula/cpp/libs/gsl/lib -I/home/leon/bebopalula/cpp/libs/gsl/include -std=c++11 ${temp_main} -Ofast -s -o ${temp_out}_${ini}.out -lgsl -lgslcblas 
    
    ./mem_usage.sh 
    screen -dmS on_${n_pop}_pop_${dir}_N_${N}_K_${K}_ini_${ini}_on_far ./${temp_out}.out $n_pop $N $K ${dir}_on 
    
done

rm $temp_main
rm $temp_globals
rm ${temp_out}.out

