#!/usr/bin/env bash 

temp_files=$(./bash_utils.sh) 
temp_globals=$(echo $temp_files | awk 'NR==1{print $1}') 
temp_main=$(echo $temp_files | awk 'NR==1{print $2}') 
temp_out=$(echo $temp_files | awk 'NR==1{print $3}') 

IF_LIF=1  
IF_BIN=0

IF_STP=1

IF_GAUSS=0
IF_RING=1

IF_GEN_CON=0
IF_SPEC=0 

IF_LOW_RANK=0
RANK=1
FIX_MAP_SEED=1

IF_HYSTERESIS=0 
HYST_J_EE=0 

sed -ie "s/ DURATION .*/ DURATION (double) 10E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_STEADY .*/ TIME_STEADY (double) 10E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_WINDOW .*/ TIME_WINDOW (double) 0.2500E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_REC .*/ TIME_REC (double) 14E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_REC_SPIKES .*/ TIME_REC_SPIKES (double) 0E3 /" "$temp_globals" ; 

sed -ie "s/ IF_LIF .*/ IF_LIF ${IF_LIF} /" "$temp_globals" ; 
sed -ie "s/ IF_BIN .*/ IF_BIN ${IF_BIN} /" "$temp_globals" ; 

sed -ie "s/ IF_STP .*/ IF_STP ${IF_STP} /" "$temp_globals" ; 

sed -ie "s/ IF_TRIALS .*/ IF_TRIALS 0 /" "$temp_globals" ; 
sed -ie "s/ IF_INI_COND .*/ IF_INI_COND 1 /" "$temp_globals" ; 

sed -ie "s/ IF_GEN_CON .*/ IF_GEN_CON ${IF_GEN_CON} /" "$temp_globals" ;
sed -ie "s/ IF_SAVE_CON_VEC .*/ IF_SAVE_CON_VEC 0 /" "$temp_globals" ;
sed -ie "s/ IF_SAVE_SPARSE_REP .*/ IF_SAVE_SPARSE_REP 0 /" "$temp_globals" ;

sed -ie "s/ IF_SPEC .*/ IF_SPEC ${IF_SPEC} /" "$temp_globals" ;
sed -ie "s/ RANK .*/ RANK ${RANK} /" "$temp_globals" ;

sed -ie "s/ IF_HYSTERESIS .*/ IF_HYSTERESIS ${IF_HYSTERESIS} /" "$temp_globals" ; 
sed -ie "s/ HYST_J_EE .*/ HYST_J_EE ${HYST_J_EE} /" "$temp_globals" ; 

# sed -ie "s/ IF_DPA .*/ IF_DPA ${IF_DPA} /" "$temp_globals" ; 
# sed -ie "s/ IF_DUAL .*/ IF_DUAL ${IF_DUAL} /" "$temp_globals" ; 
# sed -ie "s/ IF_DRT .*/ IF_DRT ${IF_DRT} /" "$temp_globals" ; 

sed -ie "s/ IF_STEP .*/ IF_STEP 1 /" "$temp_globals" ; 

sed -ie "s/ FIX_MAP_SEED .*/ FIX_MAP_SEED ${FIX_MAP_SEED} /" "$temp_globals" ; 

read n_pop N K dir kappa n_ini <<<  "$1 $2 $3 $4 $5 $6" 

sed -ie "s/ KAPPA (double) .*/ KAPPA (double) ${kappa} /" "$temp_globals" ; 
sed -ie "s/ KAPPA_1 (double) .*/ KAPPA_1 (double) ${kappa} /" "$temp_globals" ; 

for ini_cond in $(seq 1 1 $n_ini); do 

    # if (($FIX_MAP_SEED==0)); then
    # sed -ie "s/ MAP_SEED .*/ MAP_SEED ${ini_cond} /" "$temp_globals" ; 
    # fi
    
    echo "#########################################################################"
    
    ./mem_usage.sh
    ./cpu_usage.sh
    
    echo "simulation parameters:" 
    echo "n_pop ${n_pop} n_neurons ${N}0000 K ${K} ${dir} kappa ${kappa} ini_cond ${ini_cond}" 
    
    sed -ie "s/INI_COND_ID .*/INI_COND_ID ${ini_cond} /" "$temp_globals" 
    
    echo "running simulation ..." 
    echo "#########################################################################" 
    
    echo "g++ $temp_main -Ofast -s -std=c++11 -o $temp_out -lgsl -lblas" 
    g++ ${temp_main} -Ofast -s -std=c++11 -o ${temp_out} -lgsl -lblas 
    
    # ./${temp_out} $n_pop $N $K $dir 
    screen -dmS ${n_pop}_pop_${dir}_N_${N}_K_${K}_kappa_${kappa}_ini_cond_${ini_cond} ./${temp_out} $n_pop $N $K $dir 
    
done

rm $temp_files
