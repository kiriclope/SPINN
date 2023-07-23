#!/usr/bin/env bash
LANG=en_US
export LANG

read dir n_pop N K <<<  "$1 $2 $3 $4" 
read pMin dp pMax <<< "$5 ${6} ${7}"

J0=1.0
sed -ie "s/J0 .*/J0 -${J0}/" globals.h 

J1=0.0
sed -ie "s/MEAN_XI .*/MEAN_XI ${J1}/" globals.h 

for i in $(seq 1.1 0.1 2.0); do 
    
    sed -ie "s/VAR_XI .*/VAR_XI $i/" globals.h 

    (cd ../../../cuda/connectivity/; make -B)
    sleep 2
    (cd ../../../cuda/connectivity/; ./a.out)
    sleep 10
    
    echo "Compiling ..."
    g++ main.cpp -std=c++11 -Ofast -s -o loop_sigma.out
    sleep 2
    
    screen -dmS  N${N}K${K}_MEAN_XI_${J1}_VAR_XI_${i} ./loop_sigma.out $dir $n_pop $N $K
    echo N${N}K${K}_MEAN_XI_${J1}_VAR_XI_${i} "Running ..."
done
