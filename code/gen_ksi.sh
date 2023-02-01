#!/usr/bin/env bash 

read n_pop N K rank seed <<<  "$1 $2 $3 $4 $5"

sed -ie "s/ IF_LIF .*/ IF_LIF 0 /" globals.h ;
sed -ie "s/ IF_LOW_RANK .*/ IF_LOW_RANK 1 /" globals.h ;

sed -ie "s/ RANK .*/ RANK $rank /" globals.h ;

sed -ie "s/ IF_GEN_KSI .*/ IF_GEN_KSI 1 /" globals.h ;
sed -ie "s/ FIX_KSI_SEED .*/ FIX_KSI_SEED 1 /" globals.h ;
sed -ie "s/ SEED_KSI .*/ SEED_KSI ${seed} /" globals.h ;

sed -ie "s/ IF_GEN_CON .*/ IF_GEN_CON 0 /" globals.h ;
sed -ie "s/ IF_GEN_CON_FF .*/ IF_GEN_CON_FF 0 /" globals.h ;

make > /dev/null 2>&1 
./a.out $n_pop $N $K ksi
# screen -dmS gen_ksi_con ./a.out $n_pop $N $K ksi
sed -ie "s/ IF_GEN_KSI .*/ IF_GEN_KSI 0 /" globals.h ;
