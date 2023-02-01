#!/usr/bin/env bash 

temp_files=$(./bash_utils.sh)
temp_globals=$(echo $temp_files | awk 'NR==1{print $1}') 
temp_main=$(echo $temp_files | awk 'NR==1{print $2}') 
temp_out=$(echo $temp_files | awk 'NR==1{print $3}') 

IF_LIF=1 
IF_BIN=0 

IF_STP=1 
IF_GEN_CON=0 

RANK=1 
IF_SPEC=0 
FIX_MAP_SEED=1 

IF_LOW_RANK=0 
FIX_KSI_SEED=1 

sed -ie "s/ DURATION .*/ DURATION (double) 10E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_INI .*/ TIME_INI (double) 2E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_STEADY .*/ TIME_STEADY (double) 20E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_WINDOW .*/ TIME_WINDOW (double) .250E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_REC .*/ TIME_REC (double) 10E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_REC_SPIKES .*/ TIME_REC_SPIKES (double) 2E3 /" "$temp_globals" ; 

sed -ie "s/ IF_LIF .*/ IF_LIF ${IF_LIF} /" "$temp_globals" ; 
sed -ie "s/ IF_BIN .*/ IF_BIN ${IF_BIN} /" "$temp_globals" ; 

sed -ie "s/ IF_STP .*/ IF_STP ${IF_STP} /" "$temp_globals" ; 

sed -ie "s/ IF_TRIALS .*/ IF_TRIALS 1 /" "$temp_globals" ; 
sed -ie "s/ IF_INI_COND .*/ IF_INI_COND 0 /" "$temp_globals" ; 

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
sed -ie "s/ IF_DPA .*/ IF_DPA 0 /" "$temp_globals" ; 
sed -ie "s/ IF_DUAL .*/ IF_DUAL 0 /" "$temp_globals" ; 

read n_pop N K dir kappa<<<  "$1 $2 $3 $4 $5" 
read gain_min dgain gain_max <<< "$6 $7 $8" 

# ( cd ../../../cuda/connectivity/ ;
#   sed -ie "s/ n_pop .*/ n_pop ${n_pop} /" globals.h ; 
#   sed -ie "s/ N_NEURONS (unsigned long) .*/ N_NEURONS (unsigned long) ${N}0000 /" globals.h ; 
#   sed -ie "s/ K (double) .*/ K (double) ${K} /" globals.h ; 
#   sed -ie "s/ IF_SPEC .*/ IF_SPEC ${IF_SPEC} /" globals.h ; 
#   sed -ie "s/ IF_SAVE_CON_VEC .*/ IF_SAVE_CON_VEC 0 /" globals.h ; 
#   sed -ie "s/ IF_SAVE_SPARSE_REP .*/ IF_SAVE_SPARSE_REP 1 /" globals.h ; 
# ) 

sed -ie "s/ SEED_KSI (double) .*/ SEED_KSI (double) 2.0 /" "$temp_globals" 

# echo "g++ $temp_main -Ofast -s -std=c++11 -o $temp_out -lgsl -lblas" 
# g++ gen_ksi.cpp -Ofast -o ksi.out -s -std=c++11 -lgsl -lblas 
# ./ksi.out $n_pop $N $K $dir 

sed -ie "s/ KAPPA (double) .*/ KAPPA (double) ${kappa} /" "$temp_globals" ; 
sed -ie "s/ KAPPA_1 (double) .*/ KAPPA_1 (double) ${kappa} /" "$temp_globals" ; 

for trial in $(seq 0 1 20); do 
    
    sed -ie "s/ MAP_SEED .*/ MAP_SEED ${trial} /" "$temp_globals" 
    sed -ie "s/ SEED_KSI (double) .*/ SEED_KSI (double) 2.0 /" "$temp_globals" 
    
    for gain in $(seq ${gain_min} ${dgain} ${gain_max}); do 

	sed -ie "s/ GAIN (double) .*/ GAIN (double) ${gain} /" "$temp_globals" ; 
	
	echo "#########################################################################" 
	./mem_usage.sh 
	./cpu_usage.sh 
	echo "#########################################################################" 
	
	# ( cd ../../../cuda/connectivity/ ; 
	#   sed -ie "s/ KAPPA (double) .*/ KAPPA (double) ${kappa} /" globals.h ; 
	#   make > /dev/null 2>&1 ; 
	# ) 
	
	echo "simulation parameters:" 
	echo "n_pop ${n_pop} n_neurons ${N}0000 K ${K} ${dir} gain ${gain} trial ${trial}" 
	echo "#########################################################################" 
	
	# ( cd ../../../cuda/connectivity/ ; 
	#   ./a.out > /dev/null ; 
	# )
	
	sed -ie "s/ TRIAL_ID .*/ TRIAL_ID ${trial} /" "$temp_globals" 
	
	echo "g++ $temp_main -Ofast -s -std=c++11 -o $temp_out -lgsl -lblas" 
	g++ ${temp_main} -Ofast -s -std=c++11 -o ${temp_out} -lgsl -lblas
	# ./${temp_out} $n_pop $N $K $dir 	
	screen -dmS ${n_pop}_pop_${dir}_N_${N}_K_${K}_gain_${gain}_trial_${trial} ./${temp_out} $n_pop $N $K $dir 
	
	# sleep 5s 
    done
    
done

# ( cd ../../../cuda/connectivity/ ; 
#   cp globals.h.bak globals.h ; 
# )

rm $temp_files
