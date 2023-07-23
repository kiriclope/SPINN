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

sed -ie "s/ DURATION .*/ DURATION (double) 10E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_STEADY .*/ TIME_STEADY (double) 10E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_WINDOW .*/ TIME_WINDOW (double) 0.250E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_REC .*/ TIME_REC (double) 10E3 /" "$temp_globals" ; 
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
sed -ie "s/ IF_INI_COND .*/ IF_INI_COND 0 /" "$temp_globals" ; 
sed -ie "s/ IF_LOOP_M0 .*/ IF_LOOP_M0 0 /" "$temp_globals" ; 

read n_pop N K dir KAPPA <<<  "$1 $2 $3 $4 $5" 
read J0_min dJ0 J0_max <<< "$6 $7 $8" 
read Jee_min dJee Jee_max <<< "$9 ${10} ${11}" 

sed -ie "s/ KAPPA (double) .*/ KAPPA (double) ${KAPPA} /" "$temp_globals" ; 

echo "#########################################################################" 

for J0 in $(seq ${J0_min} ${dJ0} ${J0_max}); do

    sed -ie "s/ M0 .*/ M0 ${J0} /" "$temp_globals" ; 

    for Jee in $(seq ${Jee_min} ${dJee} ${Jee_max}); do 

	./mem_usage.sh 
	./cpu_usage.sh
	
	dum_J0=$(printf '%.3f' ${J0}) 
	dum_Jee=$(printf '%.2f' ${Jee}) 
	
	(cd ../parameters/${n_pop}pop/ ; 
	 echo "ext_inputs 0.5 0.25" > ${dir}_J0_${dum_J0}_Jee_${dum_Jee}.txt ; 
	 echo "J ${dum_Jee} -0.375 1.8 -0.45" >> ${dir}_J0_${dum_J0}_Jee_${dum_Jee}.txt ; 
	 echo "Tsyn 3 2 3 2" >> ${dir}_J0_${dum_J0}_Jee_${dum_Jee}.txt ; 
	)
    	
	echo "simulation parameters:" 
	echo "n_pop ${n_pop} n_neurons ${N}0000 K ${K} ${dir} kappa ${KAPPA} J0 ${dum_J0} Jee ${dum_Jee}" 
	echo "#########################################################################" 
	
	echo "g++ $temp_main -std=c++11 -Ofast -s -o $temp_out -lgsl -lblas"
	g++ ${temp_main} -std=c++11 -Ofast -s -o ${temp_out} -lgsl -lblas 
	screen -dmS ${n_pop}_pop_N_${N}_K_${K}_${dir}_kappa_${KAPPA}_J0_${dum_J0}_Jee_${dum_Jee} ./${temp_out} $n_pop $N $K ${dir}_J0_${dum_J0}_Jee_${dum_Jee} 
	
    done
done 

rm $temp_files 
