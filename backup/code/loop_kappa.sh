#!/usr/bin/env bash 

rand=$RANDOM # has to be outside

./bash_utils.sh $rand

temp_globals=$(printf 'temp_%d_globals.h' $rand) 
temp_file=$(printf 'temp_%d' $rand) 

echo "${temp_file}"

IF_LIF=1 
IF_STP=1 

IF_LOW_RANK=1 
RANK=2

FIX_KSI_SEED=1
SEED_KSI=1

sed -ie "s/ DURATION .*/ DURATION (double) 30E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_INI .*/ TIME_INI (double) 0E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_STEADY .*/ TIME_STEADY (double) 10E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_WINDOW .*/ TIME_WINDOW (double) .250E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_REC .*/ TIME_REC (double) 100E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_REC_SPIKES .*/ TIME_REC_SPIKES (double) 0E3 /" "$temp_globals" ; 

sed -ie "s/ IF_LIF .*/ IF_LIF ${IF_LIF} /" "$temp_globals" ; 

sed -ie "s/ IF_STP .*/ IF_STP ${IF_STP} /" "$temp_globals" ; 

sed -ie "s/ IF_TRIALS .*/ IF_TRIALS 1 /" "$temp_globals" ; 
sed -ie "s/ IF_INI_COND .*/ IF_INI_COND 0 /" "$temp_globals" ; 

sed -ie "s/ IF_GEN_CON .*/ IF_GEN_CON 0 /" "$temp_globals" ; 
sed -ie "s/ IF_SAVE_CON_VEC .*/ IF_SAVE_CON_VEC 0 /" "$temp_globals" ; 
sed -ie "s/ IF_SAVE_SPARSE_REP .*/ IF_SAVE_SPARSE_REP 0 /" "$temp_globals" ; 

sed -ie "s/ RANK .*/ RANK ${RANK} /" "$temp_globals" ; 

sed -ie "s/ IF_LOW_RANK .*/ IF_LOW_RANK ${IF_LOW_RANK} /" "$temp_globals" ; 
sed -ie "s/ FIX_KSI_SEED .*/ FIX_KSI_SEED ${FIX_KSI_SEED} /" "$temp_globals" ; 
sed -ie "s/ SEED_KSI .*/ SEED_KSI (float) ${SEED_KSI} /" "$temp_globals" 

sed -ie "s/ IF_STEP .*/ IF_STEP 1 /" "$temp_globals" ; 

read n_pop N K dir n_trials <<<  "$1 $2 $3 $4 $5" 
read kappa_min dkappa kappa_max <<< "$6 $7 $8" 

sed -i "s/IF_CHRISTOS .*/IF_CHRISTOS = 0.0 /" params.cfg ; 
sed -i "s/IF_DPA .*/IF_DPA = 0.0 /" params.cfg ; 
sed -i "s/IF_DUAL .*/IF_DUAL = 0.0 /" params.cfg ; 

###########################
# generate low-rank vectors
###########################
echo "generate low-rank vectors"
./gen_ksi.sh ${n_pop} ${N} ${K} ${RANK} ${SEED_KSI} 
sleep 0.01s

#####################
# compile main
#####################
make ${temp_file}.out > /dev/null 2>&1

    
for kappa in $(seq ${kappa_min} ${dkappa} ${kappa_max}); do 
    
    echo "###############################################################################" 	
    echo "# kappa ${kappa} trial ${trial} "
    
    sed -i "s/KAPPA .*/KAPPA = ${kappa} /" params.cfg ; 
    sed -i "s/KAPPA_1 .*/KAPPA_1 = ${kappa} /" params.cfg ; 
    
    echo "# generate connectivity ${n_pop} ${N} ${K} ${kappa} ${SEED_KSI}" 	
    ./gen_con.sh ${n_pop} ${N} ${K} ${kappa} ${SEED_KSI} 
    sleep 0.05s
    
    for trial in $(seq 1 1 ${n_trials}); do # trial means realisation of the Jijs
	
	sed -i "s/TRIAL_ID .*/TRIAL_ID = ${trial} /" params.cfg ; 
	sleep 0.002s
	
	echo "# run sim kappa_${kappa}_trial_${trial}_${n_pop}_pop_${dir}_N_${N}_K_${K} " 	
	screen -dmS kappa_${kappa}_trial_${trial}_${n_pop}_pop_${dir}_N_${N}_K_${K} ./${temp_file}.out $n_pop $N $K ${dir} 
	sleep 0.002s
	
    done
    
done

sed -ie "s/TRIAL_ID .*/TRIAL_ID = 1 /" params.cfg
sed -ie "s/INI_COND_ID .*/INI_COND_ID = 1 /" params.cfg 

sed -ie "s/M0 .*/M0 = 0.003 /" params.cfg ; # .3
sed -ie "s/KAPPA .*/KAPPA = 4.5 /" params.cfg ; # .3
sed -ie "s/KAPPA_1 .*/KAPPA_1 = 4.5 /" params.cfg ; # .3

rm ${temp_file}*
