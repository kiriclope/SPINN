#!/usr/bin/env bash
LANG=en_US
export LANG

read n_pop N K <<<  "$1 $2 $3" 
read pMin dp pMax <<< "$4 $5 $6"

J1=0.0
sig1=5.0

sed -ie "s/MEAN_XI .*/MEAN_XI -${J1}/" globals.h 
sed -ie "s/VAR_XI .*/VAR_XI ${sig1}/" globals.h 

sed -ie "s/IF_LEFT_RIGHT .*/IF_LEFT_RIGHT 0/" globals.h
sed -ie "s/MEAN_XI_LEFT .*/MEAN_XI_LEFT -0.0/" globals.h 
sed -ie "s/VAR_XI_LEFT .*/VAR_XI_LEFT 5.0/" globals.h 

sed -ie "s/MEAN_XI_RIGHT .*/MEAN_XI_RIGHT -0.0/" globals.h 
sed -ie "s/VAR_XI_RIGHT .*/VAR_XI_RIGHT 5.0/" globals.h 

sed -ie "s/RHO .*/RHO 1.0/" globals.h 

sed -ie "s/FIXED_XI .*/FIXED_XI 1/" globals.h 
sed -ie "s/SEED_XI .*/SEED_XI 1/" globals.h 

sed -ie "s/IF_TRIALS .*/IF_TRIALS 1/" globals.h

sed -ie "s/IF_FF .*/IF_FF 0/" globals.h
sed -ie "s/MEAN_FF .*/MEAN_FF 0.0/" globals.h 
sed -ie "s/VAR_FF .*/VAR_FF 0.0/" globals.h 

for trial in $(seq 0 1 9); do
    
    sed -ie "s/TRIAL_ID .*/TRIAL_ID ${trial}/" globals.h 

    for i in $(seq 0.1 0.1 2); do 

	echo "#########################################################################"
	echo "trial ${trial}, J0 ${i}"
	echo "#########################################################################"

	sed -ie "s/J0 .*/J0 -$i/" globals.h 

	echo "generating connectivity matrix ..."
	(cd ../../../cuda/connectivity/; make -B &>/dev/null) 
	
	(cd ../../../cuda/connectivity/; make -B &>/dev/null) 
	(cd ../../../cuda/connectivity/; ./a.out &>/dev/null) 
	
	echo "compiling ..."
	g++ main.cpp -std=c++11 -Ofast -s -o loop_J0.out
	
	screen -dmS  J0_${i}_MEAN_XI_${J1}_VAR_XI_${sig1}_trial_${trial} ./loop_J0.out $n_pop $N $K
	echo "running J0_${i}_MEAN_XI_${J1}_VAR_XI_${sig1}_trial_${trial} "
    done
    
done
