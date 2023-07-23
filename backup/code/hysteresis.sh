#!/usr/bin/env bash 

temp_files=$(./bash_utils.sh)
temp_globals=$(echo $temp_files | awk 'NR==1{print $1}')
temp_main=$(echo $temp_files | awk 'NR==1{print $2}')
temp_out=$(echo $temp_files | awk 'NR==1{print $3}')

IF_LIF=0
IF_BIN=1

IF_STP=1 

IF_GEN_CON=0 
IF_SPEC=1 
RANK=1 

IF_STIM=0 

sed -ie "s/ IF_HYSTERESIS .*/ IF_HYSTERESIS 1 /" "$temp_globals" ; 

sed -ie "s/ IF_LIF .*/ IF_LIF ${IF_LIF} /" "$temp_globals" ;
sed -ie "s/ IF_BIN .*/ IF_BIN ${IF_BIN} /" "$temp_globals" ; 

sed -ie "s/ DURATION .*/ DURATION (double) 14E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_STEADY .*/ TIME_STEADY (double) 2E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_WINDOW .*/ TIME_WINDOW (double) 2E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_REC .*/ TIME_REC (double) 14E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_REC_SPIKES .*/ TIME_REC_SPIKES (double) 0E3 /" "$temp_globals" ; 

sed -ie "s/ IF_STP .*/ IF_STP ${IF_STP} /" "$temp_globals" ; 

sed -ie "s/ IF_TRIALS .*/ IF_TRIALS 0 /" "$temp_globals" ; 
sed -ie "s/ IF_INI_COND .*/ IF_INI_COND 0 /" "$temp_globals" ; 

sed -ie "s/ IF_GEN_CON .*/ IF_GEN_CON ${IF_GEN_CON} /" "$temp_globals" ;
sed -ie "s/ IF_SAVE_CON_VEC .*/ IF_SAVE_CON_VEC 0 /" "$temp_globals" ;
sed -ie "s/ IF_SAVE_SPARSE_REP .*/ IF_SAVE_SPARSE_REP 0 /" "$temp_globals" ;

sed -ie "s/ IF_SPEC .*/ IF_SPEC ${IF_SPEC} /" "$temp_globals" ;
sed -ie "s/ RANK .*/ RANK ${RANK} /" "$temp_globals" ;

sed -ie "s/ IF_STIM .*/ IF_STIM ${IF_STIM} /" "$temp_globals" ; 

read n_pop N K dir kappa <<<  "$1 $2 $3 $4 $5" 
read hyst X_min dX X_max <<<  "$6 $7 $8 $9" 

if [ "$hyst" == "Jee" ] ; then
    HYST_JEE=1
    HYST_M0=0
fi

if [ "$hyst" == "m0" ] ; then
    HYST_JEE=0
    HYST_M0=1
fi

sed -ie "s/ HYST_X_MIN (double) .*/ HYST_X_MIN (double) ${X_min} /" "$temp_globals" ; 
sed -ie "s/ HYST_X_MAX (double) .*/ HYST_X_MAX (double) ${X_max} /" "$temp_globals" ; 
sed -ie "s/ HYST_DX (double) .*/ HYST_DX (double) ${dX} /" "$temp_globals" ; 

sed -ie "s/ KAPPA (double) .*/ KAPPA (double) ${kappa} /" "$temp_globals" ; 
sed -ie "s/ KAPPA_1 (double) .*/ KAPPA_1 (double) ${kappa} /" "$temp_globals" ; 
    
echo "#########################################################################"

./mem_usage.sh
./cpu_usage.sh

echo "simulation parameters:" 
echo "n_pop ${n_pop} n_neurons ${N}0000 K ${K} ${dir} kappa ${kappa} up" 

if (($HYST_JEE==1)) ; then
    sed -ie "s/ HYST_J_EE .*/ HYST_J_EE 1 /" "$temp_globals" ; 
    sed -ie "s/ HYST_M0 .*/ HYST_M0 0 /" "$temp_globals" ; 
fi

if (($HYST_M0==1)) ; then
    sed -ie "s/ HYST_J_EE .*/ HYST_J_EE 0 /" "$temp_globals" ; 
    sed -ie "s/ HYST_M0 .*/ HYST_M0 1 /" "$temp_globals" ; 
fi

echo "g++ $temp_main -Ofast -s -std=c++11 -o $temp_out -lgsl -lblas"
    
g++ ${temp_main} -Ofast -s -std=c++11 -o ${temp_out} -lgsl -lblas 
screen -dmS up ./${temp_out} $n_pop $N $K $dir 

echo "running simulation ..." 

echo "n_pop ${n_pop} n_neurons ${N}0000 K ${K} ${dir} kappa ${kappa} down"

if (($HYST_JEE==1)); then
    sed -ie "s/ HYST_J_EE .*/ HYST_J_EE -1 /" "$temp_globals" ; 
    sed -ie "s/ HYST_M0 .*/ HYST_M0 0 /" "$temp_globals" ; 
fi

if (($HYST_M0==1)); then
    sed -ie "s/ HYST_J_EE .*/ HYST_J_EE 0 /" "$temp_globals" ; 
    sed -ie "s/ HYST_M0 .*/ HYST_M0 -1 /" "$temp_globals" ; 
fi

echo "g++ $temp_main -Ofast -s -std=c++11 -o $temp_out -lgsl -lblas"

g++ ${temp_main} -Ofast -s -std=c++11 -o ${temp_out} -lgsl -lblas 
screen -dmS dn ./${temp_out} $n_pop $N $K $dir 

echo "running simulation ..." 
    
rm $temp_files
