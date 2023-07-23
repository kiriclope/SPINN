#!/usr/bin/env bash 

read n_pop N K kappa seed <<<  "$1 $2 $3 $4 $5"

sed -ie "s/ IF_LIF .*/ IF_LIF 0 /" globals.h ;
sed -ie "s/ IF_FF_LAYER .*/ IF_FF_LAYER 1 /" globals.h ;

dum=$(echo "print(int( 0.8 * $N *10000))" | python3) 

sed -ie "s/ N_POISSON .*/ N_POISSON 32000 /" globals.h ;

# dum=$(echo "print( 0.8 * $K )" | python3) 

sed -ie "s/ KAPPA_FF .*/ KAPPA_FF (float) ${kappa} /" globals.h ; 
sed -ie "s/ K_FF .*/ K_FF (float) $K /" globals.h ; 
sed -ie "s/ IF_POISSON_FF .*/ IF_POISSON_FF 1 /" globals.h ;

sed -ie "s/ IF_CON_DIR .*/ IF_CON_DIR 0 /" globals.h ;
sed -ie "s/ SEED_CON .*/ SEED_CON ${seed} /" globals.h ;

sed -ie "s/ IF_GEN_CON_FF .*/ IF_GEN_CON_FF 1 /" globals.h ;
sed -ie "s/ IF_SAVE_SPARSE_REP_FF .*/ IF_SAVE_SPARSE_REP_FF 1 /" globals.h ;

make 
# make > /dev/null 2>&1
./a.out $n_pop $N $K con_ff
# screen -dmS gen_ff_con ./a.out $n_pop $N $K con_ff

sed -ie "s/ IF_CON_DIR .*/ IF_CON_DIR 0 /" globals.h ;
sed -ie "s/ IF_GEN_CON_FF .*/ IF_GEN_CON_FF 0 /" globals.h ;
sed -ie "s/ IF_SAVE_SPARSE_REP_FF .*/ IF_SAVE_SPARSE_REP_FF 0 /" globals.h ;
