#!/usr/bin/env bash 

temp_files=$(./bash_utils.sh)
temp_globals=$(echo $temp_files | awk 'NR==1{print $1}')
temp_main=$(echo $temp_files | awk 'NR==1{print $2}')
temp_out=$(echo $temp_files | awk 'NR==1{print $3}')

IF_LIF=1 
IF_BIN=0

IF_STP=0

IF_GEN_CON=1 
IF_SPEC=0

sed -ie "s/ DURATION .*/ DURATION (double) 20E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_STEADY .*/ TIME_STEADY (double) 2E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_WINDOW .*/ TIME_WINDOW (double) 1E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_REC .*/ TIME_REC (double) 15E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_REC_SPIKES .*/ TIME_REC_SPIKES (double) 2E3 /" "$temp_globals" ; 

sed -ie "s/ IF_LIF .*/ IF_LIF ${IF_LIF} /" "$temp_globals" ; 
sed -ie "s/ IF_BIN .*/ IF_BIN ${IF_BIN} /" "$temp_globals" ; 

sed -ie "s/ IF_GEN_CON .*/ IF_GEN_CON ${IF_GEN_CON} /" "$temp_globals" ; 
sed -ie "s/ IF_SAVE_CON_VEC .*/ IF_SAVE_CON_VEC 0 /" "$temp_globals" ; 
sed -ie "s/ IF_SAVE_SPARSE_REP .*/ IF_SAVE_SPARSE_REP 0 /" "$temp_globals" ; 

sed -ie "s/ IF_SPEC .*/ IF_SPEC ${IF_SPEC} /" "$temp_globals" ; 

sed -ie "s/ IF_STP .*/ IF_STP ${IF_STP} /" "$temp_globals" ; 

sed -ie "s/ IF_TRIALS .*/ IF_TRIALS 0 /" "$temp_globals" ; 
sed -ie "s/ IF_INI_COND .*/ IF_INI_COND 0 /" "$temp_globals" ; 
sed -ie "s/ IF_HYSTERESIS .*/ IF_HYSTERESIS 0 /" "$temp_globals" ; 

read n_pop N dir kappa <<<  "$1 $2 $3 $4" 
read K_min dK K_max <<< "$5 $6 $7"

sed -ie "s/ KAPPA (double) .*/ KAPPA (double) ${kappa} /" "$temp_globals" ; 

for K in $(seq ${K_min} ${dK} ${K_max}); do 
    
    ./mem_usage.sh 
    ./cpu_usage.sh
    
    echo "n_pop ${n_pop} N $N K $K kappa $kappa" 
    echo "#########################################################################"
    
    # echo "generating connectivity matrix ..."    

    # (cd ../../../cuda/connectivity/ ; 
    #  sed -ie "s/ K .*/ K ${K}.0/" globals.h ; 
    #  make -B &>/dev/null ; 
    #  ./a.out &>/dev/null 
    # ) 

    echo "g++ $temp_main -std=c++11 -Ofast -s -o $temp_out -lgsl -lblas"
    g++ ${temp_main} -std=c++11 -Ofast -s -o ${temp_out} -lgsl -lblas 
    screen -dmS  ${n_pop}_pop_${dir}_N_${N}_K_$K ./${temp_out} $n_pop $N $K $dir 
        
done
