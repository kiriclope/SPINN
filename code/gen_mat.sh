#!/usr/bin/env bash
read n_pop N K pop_size<<<  "$1 $2 $3 $4" 

(cd ../../../cuda/connectivity/; sed -ie "s/n_pop .*/n_pop ${n_pop}/" globals.h) 
(cd ../../../cuda/connectivity/; sed -ie "s/N_NEURONS .*/N_NEURONS ${N}0000ULL/" globals.h) 
(cd ../../../cuda/connectivity/; sed -ie "s/ K .*/ K ${K}.0/" globals.h) 
(cd ../../../cuda/connectivity/; sed -ie "s/popSize .*/popSize ${pop_size}/" globals.h) 

echo "generating connectivity matrix ..."
(cd ../../../cuda/connectivity/; make -B &>/dev/null) 
(cd ../../../cuda/connectivity/; ./a.out &>/dev/null) 
