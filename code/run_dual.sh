#!/usr/bin/env bash 

read n_pop N K dir kappa seed <<< "$1 $2 $3 $4 $5 $6" 

# ./gen_ksi.sh ${n_pop} ${N} ${K} 2 ${seed}
# sleep 5s

# ./gen_con.sh ${n_pop} ${N} ${K} ${kappa} ${seed}
# sleep 5s

sed -i "s/ IF_LIF .*/ IF_LIF 1 /" globals.h ; 
sed -i "s/ SEED_KSI .*/ SEED_KSI ${seed} /" globals.h ; 

make 

sed -i "s/KAPPA .*/KAPPA = ${kappa} /" params.cfg ; 
sed -i "s/KAPPA_1 .*/KAPPA_1 = ${kappa} /" params.cfg ; 

sed -i "s/IF_CHRISTOS .*/IF_CHRISTOS = 0 /" params.cfg ;
sed -i "s/IF_DPA .*/IF_DPA = 0 /" params.cfg ;
sed -i "s/IF_DUAL .*/IF_DUAL = 1 /" params.cfg ;

sed -i "s/SAMPLE .*/SAMPLE = 0 /" params.cfg ;
sleep 0.001s 

screen -dmS a_${n_pop}_pop_${dir}_N_${N}_K_${K}_dpa ./a.out $n_pop $N $K ${dir}
sleep 0.001s 

# sed -i "s/SAMPLE .*/SAMPLE = 1 /" params.cfg ;

./a.out $n_pop $N $K ${dir}
# screen -dmS b_${n_pop}_pop_${dir}_N_${N}_K_${K}_dpa ./a.out $n_pop $N $K ${dir}
# sleep 0.001s 

# sed -i "s/IF_DPA .*/IF_DPA = 0.0 /" params.cfg ; 
# sed -i "s/IF_DUAL .*/IF_DUAL = 1.0 /" params.cfg ;

# sed -i "s/SAMPLE .*/SAMPLE = 0 /" params.cfg ;

# screen -dmS a_${n_pop}_pop_${dir}_N_${N}_K_${K}_dual ./a.out $n_pop $N $K ${dir}
# sleep 0.001s 

# sed -i "s/SAMPLE .*/SAMPLE = 1 /" params.cfg ;

# screen -dmS b_${n_pop}_pop_${dir}_N_${N}_K_${K}_dual ./a.out $n_pop $N $K ${dir}
# sleep 0.001s 

make clean
