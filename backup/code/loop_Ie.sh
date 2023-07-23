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

IF_STIM=0 

sed -ie "s/ IF_LIF .*/ IF_LIF ${IF_LIF} /" "$temp_globals" ; 
sed -ie "s/ IF_BIN .*/ IF_BIN ${IF_BIN} /" "$temp_globals" ; 

sed -ie "s/ DURATION .*/ DURATION (double) 2E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_STEADY .*/ TIME_STEADY (double) 10E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_WINDOW .*/ TIME_WINDOW (double) 2E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_REC .*/ TIME_REC (double) 2E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_REC_SPIKES .*/ TIME_REC_SPIKES (double) 0E3 /" "$temp_globals" ; 

sed -ie "s/ IF_STP .*/ IF_STP ${IF_STP} /" "$temp_globals" ; 

sed -ie "s/ IF_GEN_CON .*/ IF_GEN_CON ${IF_GEN_CON} /" "$temp_globals" ; 
sed -ie "s/ IF_SAVE_CON_VEC .*/ IF_SAVE_CON_VEC 0 /" "$temp_globals" ;
sed -ie "s/ IF_SAVE_SPARSE_REP .*/ IF_SAVE_SPARSE_REP 0 /" "$temp_globals" ;

sed -ie "s/ IF_SPEC .*/ IF_SPEC ${IF_SPEC} /" "$temp_globals" ; 
sed -ie "s/ RANK .*/ RANK ${RANK} /" "$temp_globals" ; 

sed -ie "s/ IF_STIM .*/ IF_STIM ${IF_STIM} /" "$temp_globals" ; 
sed -ie "s/ IF_HYSTERESIS .*/ IF_HYSTERESIS 0 /" "$temp_globals" ;

sed -ie "s/ IF_TRIALS .*/ IF_TRIALS 0 /" "$temp_globals" ; 
sed -ie "s/ IF_INI_COND .*/ IF_INI_COND 1 /" "$temp_globals" ; 

read n_pop N K dir KAPPA <<<  "$1 $2 $3 $4 $5" 
read Ie_min dIe Ie_max <<< "$6 $7 $8" 

sed -ie "s/ KAPPA (double) .*/ KAPPA (double) ${KAPPA} /" "$temp_globals" ; 

echo "#########################################################################" 

for ini in $(seq 0 1 10); do 
    
    sed -ie "s/INI_COND_ID .*/INI_COND_ID ${ini} /" "$temp_globals"  
    
    for Ie in $(seq ${Ie_min} ${dIe} ${Ie_max}); do 

	./mem_usage.sh 
	./cpu_usage.sh
	
	dum=$(printf '%.2f' ${Ie}) 
	
	(cd ../parameters/${n_pop}pop/ ; 
	 echo "ext_inputs ${dum} 0.25" > ${dir}_Ie_${dum}.txt ; 
	 echo "J 2.00 -0.75 1.8 -0.9" >> ${dir}_Ie_${dum}.txt ; 
	 echo "Tsyn 3 2 3 2" >> ${dir}_Ie_${dum}.txt ; 
	)
    	
	echo "simulation parameters:" 
	echo "n_pop ${n_pop} n_neurons ${N}0000 K ${K} ${dir} kappa ${KAPPA} Ie ${dum} ini ${ini}" 
	echo "#########################################################################" 
	
	echo "g++ $temp_main -std=c++11 -Ofast -s -o $temp_out -lgsl -lblas"
	g++ ${temp_main} -std=c++11 -Ofast -s -o ${temp_out} -lgsl -lblas 
	screen -dmS  ${n_pop}_pop_N_${N}_K_${K}_${dir}_kappa_${KAPPA}_Ie_${dum}_ini_${ini} ./${temp_out} $n_pop $N $K ${dir}_Ie_${dum} 
	
    done
done 

rm $temp_files 
