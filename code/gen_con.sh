#!/usr/bin/env bash 

read n_pop N K kappa seed <<<  "$1 $2 $3 $4 $5" 

# generating connectivity with cuda 
( cd ../../../cuda/connectivity/ ; 
  sed -ie "s/ n_pop .*/ n_pop ${n_pop} /" globals.h ;
  sed -ie "s/ N_NEURONS .*/ N_NEURONS (unsigned long) ${N}0000 /" globals.h ;
  sed -ie "s/ K .*/ K (float) ${K} /" globals.h ;
  
  sed -ie "s/ IF_RING .*/ IF_RING 1 /" globals.h ;
  
  sed -ie "s/ IF_LOW_RANK .*/ IF_LOW_RANK 0 /" globals.h ; 
  sed -ie "s/ RANK .*/ RANK 2 /" globals.h ; 
  sed -ie "s/ SEED_KSI .*/ SEED_KSI $seed /" globals.h ; 
  
  sed -ie "s/ KAPPA (float) .*/ KAPPA (float) ${kappa} /" globals.h ; 
  sed -ie "s/ IF_SAVE_SPARSE_REP .*/ IF_SAVE_SPARSE_REP 1 /" globals.h ; 
  sed -ie "s/ IF_CON_DIR .*/ IF_CON_DIR 0 /" globals.h ; 
  sed -ie "s/ SEED_CON .*/ SEED_CON (float) 3.0 /" globals.h ; 
  sed -ie "s/ IF_REMAP_SQRT_K .*/ IF_REMAP_SQRT_K 0 /" globals.h ;
  
  make > /dev/null 2>&1 ; 
  ./a.out 
) 
