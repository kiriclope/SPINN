#!/usr/bin/env bash 

rand=$RANDOM # has to be outside

./bash_utils.sh $rand

temp_globals=$(printf 'temp_%d_globals.h' $rand) 
temp_file=$(printf 'temp_%d' $rand) 

IF_LIF=1 
IF_BIN=0 

IF_STP=1
IF_GEN_CON=0
IF_CUDA=0

RANK=1 
IF_GAUSS=0 
IF_RING=1 

IF_SPEC=0
FIX_MAP_SEED=1 

IF_LOW_RANK=0 
FIX_KSI_SEED=1 

sed -ie "s/ DURATION .*/ DURATION (float) 8E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_INI .*/ TIME_INI (float) 0E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_STEADY .*/ TIME_STEADY (float) 6E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_WINDOW .*/ TIME_WINDOW (float) .050E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_REC .*/ TIME_REC (float) 60E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_REC_SPIKES .*/ TIME_REC_SPIKES (float) 0E3 /" "$temp_globals" ; 

sed -ie "s/ IF_LIF .*/ IF_LIF ${IF_LIF} /" "$temp_globals" ; 
sed -ie "s/ IF_BIN .*/ IF_BIN ${IF_BIN} /" "$temp_globals" ; 

sed -ie "s/ IF_STP .*/ IF_STP ${IF_STP} /" "$temp_globals" ; 

sed -ie "s/ IF_TRIALS .*/ IF_TRIALS 1 /" "$temp_globals" ; 
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

sed -ie "s/ IF_STEP .*/ IF_STEP 0 /" "$temp_globals" ; 

read n_pop N K dir n_trials n_ini <<<  "$1 $2 $3 $4 $5 $6" 

make ${temp_file}.out

for ini in $(seq 1 1 $n_ini); do 
    
    # sed -ie "s/ INI_COND_ID .*/ INI_COND_ID ${ini} /" "$temp_globals"     
    sed -ie "s/INI_COND_ID .*/INI_COND_ID = ${ini} /" params.cfg 
    
    if (($IF_CUDA==1)); then 
	sed -ie "s/ IF_CON_DIR .*/ IF_CON_DIR 1 /" "$temp_globals" ; 
	# sed -ie "s/ SEED_CON .*/ SEED_CON (float) ${ini} /" "$temp_globals" ; 
	sed -ie "s/ SEED_CON .*/ SEED_CON (float) 3 /" "$temp_globals" ; 	
    fi
    
    for trial in $(seq 1 1 $n_trials); do
	
	echo "#########################################################################"
	
	dum=$(echo "print(${trial} / ${n_trials})" | python3) 
	dum2=$(echo "print(0.25 + ${trial} / ${n_trials})" | python3) 
	
	echo "simulation parameters:" 
	echo "n_pop ${n_pop} n_neurons ${N}0000 K ${K} ${dir} trial ${trial} phi ${dum} ini ${ini}" 
	echo "#########################################################################" 
	
	sed -ie "s/PHI_CUE .*/PHI_CUE = $dum /" params.cfg ; 
	sed -ie "s/PHI_ERASE .*/PHI_ERASE = $dum2 /" params.cfg ; 
	
	sed -ie "s/TRIAL_ID .*/TRIAL_ID = ${trial} /" params.cfg
	
	./mem_usage.sh 
	./cpu_usage.sh 
	
	sed -ie "s/M0 .*/M0 = 0.0027 /" params.cfg ; # .3
	sed -i "s/VlE .*/VlE = -0.5 /" params.cfg  ; 
	# sed -ie "s/NU_FF .*/NU_FF = 0.0 /" params.cfg ; # .3 
	sleep 0.002s
	# sed -ie "s/ IF_CON_DIR .*/ IF_CON_DIR 0 /" "$temp_globals" ; 
	# sed -ie "s/ SEED_CON .*/ SEED_CON (float) 3 /" "$temp_globals" ; 
	
	screen -dmS off_ini_${ini}_trial_${trial}_${n_pop}_pop_${dir}_N_${N}_K_${K} ./${temp_file}.out $n_pop $N $K ${dir}_off 
	sleep 0.002s 

	./mem_usage.sh 
	./cpu_usage.sh 
	
	# sed -ie "s/M0 .*/M0 = 0.004 /" params.cfg ; # .3 
	# sed -ie "s/NU_FF .*/NU_FF = 0.00125 /" params.cfg ; # .3 
	sed -i "s/M0 .*/M0 = 0.0033 /" params.cfg ; 
	sed -i "s/VlE .*/VlE = 1.0 /" params.cfg  ; 
	sleep 0.002s 
	
	screen -dmS on_ini_${ini}_trial_${trial}_${n_pop}_pop_${dir}_N_${N}_K_${K} ./${temp_file}.out $n_pop $N $K ${dir}_on 
	sleep 0.002s 
	
	# g++ ${temp_main} -Ofast -s -std=c++11 -o ${temp_out} -lgsl -lblas 
	# ./${temp_out} $n_pop $N $K $dir 
	# screen -dmS ${n_pop}_pop_${dir}_N_${N}_K_${K}_trial_${trial}_ini_${ini} ./${temp_out}_${trial}_${ini}.out $n_pop $N $K $dir 
	# screen -dmS ${n_pop}_pop_${dir}_N_${N}_K_${K}_trial_${trial}_ini_${ini}_off srun --priority="TOP" ./${temp_out}_${ini}_${trial}.out $n_pop $N $K albert_off 
	# screen -dmS ${n_pop}_pop_${dir}_N_${N}_K_${K}_trial_${trial}_ini_${ini}_on srun --priority="TOP" ./${temp_out}_${ini}_${trial}.out $n_pop $N $K albert_on
	
    done
        
done

sleep 5s

sed -ie "s/PHI_CUE .*/PHI_CUE = 0.25 /" params.cfg 
sed -ie "s/PHI_ERASE .*/PHI_ERASE = 0.75 /" params.cfg 

sed -ie "s/TRIAL_ID .*/TRIAL_ID = 1 /" params.cfg
sed -ie "s/INI_COND_ID .*/INI_COND_ID = 1 /" params.cfg 

sed -ie "s/M0 .*/M0 = 0.003 /" params.cfg ; # .3

rm ${temp_file}*
