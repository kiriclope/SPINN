#!/usr/bin/env bash
LANG=en_US
export LANG

read n_pop N K <<<  "$1 $2 $3" 
read pMin dp pMax <<< "$4 $5 $6"

J0=0.5 
J1=0.0 
sig1=5.0 

sed -ie "s/J0 .*/J0 -${J0}/" globals.h 

sed -ie "s/MEAN_XI .*/MEAN_XI -${J1}/" globals.h 
sed -ie "s/VAR_XI .*/VAR_XI ${sig1}/" globals.h 

sed -ie "s/FIXED_XI .*/FIXED_XI 1/" globals.h 
sed -ie "s/SEED_XI .*/SEED_XI 1/" globals.h 

sed -ie "s/IF_TRIALS .*/IF_TRIALS 1/" globals.h 

sed -ie "s/IF_FF .*/IF_FF 1/" globals.h 
sed -ie "s/FIXED_FF .*/FIXED_FF 1/" globals.h 

for trial in $(seq 0 1 3); do 
    
    sed -ie "s/TRIAL_ID .*/TRIAL_ID ${trial}/" globals.h 
    
    for var in $(seq 0.0 1.0 10.0); do 

	echo "trial ${trial}, ff_var_ortho ${var}"
	echo "#########################################################################"

	sed -ie "s/VAR_ORTHO .*/VAR_ORTHO ${var}/" globals.h 

	echo "generating connectivity matrix ..."
	(cd ../../../cuda/connectivity/; make -B &>/dev/null) 
	(cd ../../../cuda/connectivity/; ./a.out &>/dev/null ) 
	
	echo "compiling ..."
	g++ main.cpp -std=c++11 -Ofast -s -o loop_ff_ortho.out	

	screen -dmS  J0_${J0}_MEAN_XI_${J1}_VAR_XI_${sig1}_trial_${trial}_ortho_${var} ./loop_ff_ortho.out $dir $n_pop $N $K
	echo "running J0_${J0}_MEAN_XI_${J1}_VAR_XI_${sig1}_trial_${trial}_ortho_${var}"
	echo "#########################################################################"
    done

done
