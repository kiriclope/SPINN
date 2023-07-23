#!/usr/bin/env bash 

read n_pop N K dir kappa <<< "$1 $2 $3 $4 $5" 
seed=1

# echo "generate FF connectiviy"
# ./gen_con_ff.sh ${n_pop} ${N} ${K} .5 ${seed}
# sleep 5s

# echo "generate connectiviy"
# ./gen_con.sh ${n_pop} ${N} ${K} ${kappa} ${seed}
# sleep 5s

sed -i "s/ IF_LIF .*/ IF_LIF 1 /" globals.h ; 

echo "compile"
make 
sed -i "s/M0 .*/M0 = 0.0027 /" params.cfg  ; 

sed -i "s/NU_FF .*/NU_FF = 0.0 /" params.cfg ;
sed -i "s/KAPPA .*/KAPPA = ${kappa} /" params.cfg ; 
# sed -i "s/VlE .*/VlE = -0.5 /" params.cfg  ; 
sleep 0.001s 

echo "run sim off"
# ./a.out $n_pop $N $K ${dir}_off
screen -dmS a_${n_pop}_pop_${dir}_N_${N}_K_${K}_off ./a.out $n_pop $N $K ${dir}_off
sleep 0.001s 

# sed -i "s/M0 .*/M0 = 0.0033 /" params.cfg ; 
# sed -i "s/VlE .*/VlE = 1.0 /" params.cfg  ; 
# sleep 0.001s 

echo "run sim on"
screen -dmS b_${n_pop}_pop_${dir}_N_${N}_K_${K}_on ./a.out $n_pop $N $K ${dir}_on 
sleep 0.001s 

# sed -i "s/VlE .*/VlE = 0.0 /" params.cfg  ;

make clean
