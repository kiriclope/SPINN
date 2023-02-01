#!/usr/bin/env bash 

temp_files=$(./bash_utils.sh)
temp_globals=$(echo $temp_files | awk 'NR==1{print $1}')
temp_main=$(echo $temp_files | awk 'NR==1{print $2}')
temp_out=$(echo $temp_files | awk 'NR==1{print $3}')

IF_LIF=1 
IF_BIN=0 

IF_STP=1 

IF_GEN_CON=0 
IF_SPEC=1 
RANK=1 
FIX_MAP_SEED=0 

sed -ie "s/ IF_LOOP_M0 .*/ IF_LOOP_M0 1 /" "$temp_globals" ; 

sed -ie "s/ DURATION .*/ DURATION (double) 10E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_STEADY .*/ TIME_STEADY (double) 10E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_WINDOW .*/ TIME_WINDOW (double) 2E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_REC .*/ TIME_REC (double) 10E3 /" "$temp_globals" ; 
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
sed -ie "s/ FIX_MAP_SEED .*/ FIX_MAP_SEED ${FIX_MAP_SEED} /" "$temp_globals" ; 

sed -ie "s/ IF_STIM .*/ IF_STIM 0 /" "$temp_globals" ; 
sed -ie "s/ IF_HYSTERESIS .*/ IF_HYSTERESIS 0 /" "$temp_globals" ; 

read n_pop N K dir kappa <<<  "$1 $2 $3 $4 $5" 
read m0_min dm0 m0_max <<< "$6 $7 $8"

sed -ie "s/ KAPPA (double) .*/ KAPPA (double) ${kappa} /" "$temp_globals" ; 
sed -ie "s/ KAPPA_1 (double) .*/ KAPPA_1 (double) ${kappa} /" "$temp_globals" ; 

for ini in $(seq 0 1 10); do 
    
    for m0 in $(seq ${m0_min} ${dm0} ${m0_max}); do 
	
	echo "#########################################################################" 
	./mem_usage.sh 
	./cpu_usage.sh 
	echo "#########################################################################" 
	
	sed -ie "s/ M0 (double) .*/ M0 (double) ${m0}E-3 /" "$temp_globals" ; 
	
	echo "simulation parameters:" 
	echo "n_pop ${n_pop} n_neurons ${N}0000 K ${K} ${dir} m0 ${m0} ini ${ini}" 
	echo "#########################################################################" 
		
	sed -ie "s/INI_ID .*/INI_ID ${ini} /" "$temp_globals" 
	
	echo "g++ $temp_main -Ofast -s -std=c++11 -o $temp_out -lgsl -lblas"
	g++ ${temp_main} -Ofast -s -std=c++11 -o ${temp_out} -lgsl -lblas 
	screen -dmS ${n_pop}_pop_${dir}_N_${N}_K_${K}_m0_${m0}_ini_${ini} ./${temp_out} $n_pop $N $K $dir 
	
    done
    
done

rm $temp_files
