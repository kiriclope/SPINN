#!/bin/bash

read rand <<< "$1" 

temp_globals=$(printf 'temp_%d_globals.h' $rand) ;
temp_file=$(printf 'temp_%d.cpp' $rand) ;

cp main.cpp $temp_file ;
cp globals.h $temp_globals ;

sed -i 's/ "globals.h" .*/ "'"$temp_globals"'" /' "$temp_file" ; 

