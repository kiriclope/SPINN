#!/usr/bin/env bash 

temp_files=$(./bash_utils.sh)
temp_globals=$(echo $temp_files | awk 'NR==1{print $1}') 
temp_main=$(echo $temp_files | awk 'NR==1{print $2}') 
temp_out=$(echo $temp_files | awk 'NR==1{print $3}') 

IF_LIF=1
IF_BIN=0 

IF_STP=1 
if (($IF_STP==1)); then 
    sed -ie "s/ M0 (double) .*/ M0 (double) 0.05 /" "$temp_globals" ; 
fi
if (($IF_STP==0)); then 
    sed -ie "s/ M0 (double) .*/ M0 (double) 0.01 /" "$temp_globals" ; 
fi

if (($IF_LIF==1)); then 
    if (($IF_STP==1)); then 
	sed -ie "s/ M0 (double) .*/ M0 (double) 0.0025 /" "$temp_globals" ; 
    fi
    if (($IF_STP==0)); then 
	sed -ie "s/ M0 (double) .*/ M0 (double) 0.005 /" "$temp_globals" ; 
    fi
fi

IF_GEN_CON=0 

RANK=1 
IF_SPEC=0 
FIX_MAP_SEED=1 

IF_LOW_RANK=1 
FIX_KSI_SEED=1 

sed -ie "s/ DURATION .*/ DURATION (double) 10E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_STEADY .*/ TIME_STEADY (double) 30E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_WINDOW .*/ TIME_WINDOW (double) .250E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_REC .*/ TIME_REC (double) 10E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_REC_SPIKES .*/ TIME_REC_SPIKES (double) 2E3 /" "$temp_globals" ; 

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

sed -ie "s/ IF_STIM .*/ IF_STIM 0 /" "$temp_globals" ; 
sed -ie "s/ IF_HYSTERESIS .*/ IF_HYSTERESIS 0 /" "$temp_globals" ; 

read n_pop N K dir <<<  "$1 $2 $3 $4" 
read ksi_var_min dksi_var ksi_var_max <<< "$5 $6 $7" 

for ini_cond in $(seq 0 1 20); do 
    
    sed -ie "s/ SEED_KSI (double) .*/ SEED_KSI (double) 2.0 /" "$temp_globals" 
    
    for kappa in $(seq ${kappa_min} ${dkappa} ${kappa_max}); do 
	
	echo "#########################################################################" 
	./mem_usage.sh 
	./cpu_usage.sh 
	echo "#########################################################################" 
	
	sed -ie "s/ KAPPA (double) .*/ KAPPA (double) ${kappa} /" "$temp_globals" ; 
	sed -ie "s/ KAPPA_1 (double) .*/ KAPPA_1 (double) ${kappa} /" "$temp_globals" ; 
	
	echo "simulation parameters:" 
	echo "n_pop ${n_pop} n_neurons ${N}0000 K ${K} ${dir} kappa ${kappa} ini_cond ${ini_cond}" 
	echo "#########################################################################" 
	
	sed -ie "s/ INI_COND_ID .*/ INI_COND_ID ${ini_cond} /" "$temp_globals" 
	
	echo "g++ $temp_main -Ofast -s -std=c++11 -o $temp_out -lgsl -lblas" 
	g++ ${temp_main} -Ofast -s -std=c++11 -o ${temp_out} -lgsl -lblas 
	screen -dmS ${n_pop}_pop_${dir}_N_${N}_K_${K}_kappa_${kappa}_ini_cond_${ini_cond} ./${temp_out} $n_pop $N $K $dir 
	
	# sleep 5s 
    done
    
done

# ( cd ../../../cuda/connectivity/ ; 
#   cp globals.h.bak globals.h ; 
# )

rm $temp_files
