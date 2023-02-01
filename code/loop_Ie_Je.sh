#!/usr/bin/env bash 
LANG=en_US
export LANG

IF_LIF=1 
IF_BIN=0 

IF_GEN_CON=1
IF_STP=1  
IF_SPEC=1 
IF_RANK_2=0 
IF_STIM=0 

read n_pop N dir <<<  "$1 $2 $3" 
read K_min dK K_max <<< "$4 $5 $6" 

for K in $(seq ${K_min} ${dK} ${K_max}); do 
    
    echo "#########################################################################"
    echo "n_pop ${n_pop} N $N K $K" 
    echo "#########################################################################"

    # echo "generating connectivity matrix ..."    

    # (cd ../../../cuda/connectivity/ ; 
    #  sed -ie "s/ K .*/ K ${K}.0/" globals.h ; 
    #  make -B &>/dev/null ; 
    #  ./a.out &>/dev/null 
    # ) 
    
    echo "ext_inputs ${Ie} 0.25" > Ie_${}_Jee_
    echo "compiling ..." 
    g++ lif.cpp -std=c++11 -Ofast -s -o loop_K.out 
    screen -dmS  ${n_pop}_pop_${dir}_N_${N}_K_$K ./loop_K.out $n_pop $N $K $dir 
    
    echo "running ${n_pop}_pop_${dir}_N_${N}_K_$K" 
    
done 
