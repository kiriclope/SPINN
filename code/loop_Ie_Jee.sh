#!/usr/bin/env bash

temp_files=$(./bash_utils.sh)
temp_globals=$(echo $temp_files | awk 'NR==1{print $1}')
temp_main=$(echo $temp_files | awk 'NR==1{print $2}')
temp_out=$(echo $temp_files | awk 'NR==1{print $3}')

IF_LIF=1 
IF_BIN=0 

IF_STP=1 

IF_GEN_CON=0 
IF_SPEC=0
IF_RING=1
IF_GAUSS=0
RANK=1 

sed -ie "s/ IF_LIF .*/ IF_LIF ${IF_LIF} /" "$temp_globals" ; 
sed -ie "s/ IF_BIN .*/ IF_BIN ${IF_BIN} /" "$temp_globals" ; 

sed -ie "s/ DURATION .*/ DURATION (double) 30E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_STEADY .*/ TIME_STEADY (double) 30E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_WINDOW .*/ TIME_WINDOW (double) .5E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_REC .*/ TIME_REC (double) 30E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_REC_SPIKES .*/ TIME_REC_SPIKES (double) 0E3 /" "$temp_globals" ; 

sed -ie "s/ IF_STP .*/ IF_STP ${IF_STP} /" "$temp_globals" ; 

sed -ie "s/ IF_GEN_CON .*/ IF_GEN_CON ${IF_GEN_CON} /" "$temp_globals" ; 
sed -ie "s/ IF_SAVE_CON_VEC .*/ IF_SAVE_CON_VEC 0 /" "$temp_globals" ;
sed -ie "s/ IF_SAVE_SPARSE_REP .*/ IF_SAVE_SPARSE_REP 0 /" "$temp_globals" ;

sed -ie "s/ IF_SPEC .*/ IF_SPEC ${IF_SPEC} /" "$temp_globals" ; 
sed -ie "s/ RANK .*/ RANK ${RANK} /" "$temp_globals" ; 

sed -ie "s/ IF_HYSTERESIS .*/ IF_HYSTERESIS 0 /" "$temp_globals" ;

sed -ie "s/ IF_TRIALS .*/ IF_TRIALS 0 /" "$temp_globals" ; 
sed -ie "s/ IF_INI_COND .*/ IF_INI_COND 0 /" "$temp_globals" ; 

read n_pop N K dir KAPPA <<<  "$1 $2 $3 $4 $5" 
read Ie_min dIe Ie_max <<< "$6 $7 $8" 
read Jee_min dJee Jee_max <<< "$9 ${10} ${11}" 

sed -ie "s/ KAPPA (double) .*/ KAPPA (double) ${KAPPA} /" "$temp_globals" ; 

echo "#########################################################################" 

for Ie in $(seq ${Ie_min} ${dIe} ${Ie_max}); do         
    for Jee in $(seq ${Jee_min} ${dJee} ${Jee_max}); do 

	./mem_usage.sh 
	./cpu_usage.sh
	
	dum_Ie=$(printf '%.2f' ${Ie}) 
	dum_Jee=$(printf '%.2f' ${Jee}) 
	
	(cd ../parameters/${n_pop}pop/ ; 
	 echo "ext_inputs ${dum_Ie} 540" > ${dir}_Ie_${dum_Ie}_Jee_${dum_Jee}.txt ; 
	 echo "J ${dum_Jee} -140 120 -90" >> ${dir}_Ie_${dum_Ie}_Jee_${dum_Jee}.txt ; 
	 echo "Tsyn 3 2 3 2" >> ${dir}_Ie_${dum_Ie}_Jee_${dum_Jee}.txt ; 
	)
    	
	echo "simulation parameters:" 
	echo "n_pop ${n_pop} n_neurons ${N}0000 K ${K} ${dir} kappa ${KAPPA} Ie ${dum_Ie} Jee ${dum_Jee}" 
	echo "#########################################################################" 
	
	echo "g++ $temp_main -std=c++11 -Ofast -s -o $temp_out -lgsl -lblas"
	g++ ${temp_main} -std=c++11 -Ofast -s -o ${temp_out} -lgsl -lblas 
	screen -dmS ${n_pop}_pop_N_${N}_K_${K}_${dir}_kappa_${KAPPA}_Ie_${dum_Ie}_Jee_${dum_Jee} ./${temp_out} $n_pop $N $K ${dir}_Ie_${dum_Ie}_Jee_${dum_Jee} 
	
    done
done 

rm $temp_files 
