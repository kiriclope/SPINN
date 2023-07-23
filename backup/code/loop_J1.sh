#!/usr/bin/env bash
LANG=en_US
export LANG

read dir n_pop N K <<<  "$1 $2 $3 $4" 
read pMin dp pMax <<< "$5 ${6} ${7}"

J0=1.0
sed -ie "s/J0 .*/J0 -${J0}/" globals.h 

sig1=0.5
sed -ie "s/VAR_XI .*/VAR_XI ${sig1}/" globals.h 

for i in $(seq 0.0 0.1 1.0); do 
    
    sed -ie "s/MEAN_XI .*/MEAN_XI $i/" globals.h 
    echo "MEAN_XI_$i"
    
    (cd ../../../cuda/connectivity/; make -B)
    sleep 2
    (cd ../../../cuda/connectivity/; ./a.out)
    sleep 10

    echo "Compiling ..."
    g++ main.cpp -std=c++11 -Ofast -s -o loop_J1.out
    sleep 2
    
    screen -dmS  N${N}K${K}_MEAN_XI_${i} ./loop_J1.out $dir $n_pop $N $K
    echo N${N}K${K}_MEAN_XI_${i} "Running ..."
done
