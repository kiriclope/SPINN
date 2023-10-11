#!/usr/bin/env bash 

calc(){ awk "BEGIN { print "$*" *100 }"; } 

cpu_use () {

    n_cpu=$(nproc --all) ;
    n_threads=$(top -bn1 | awk 'NR==2{print $4}') ;
    
    cpu_usage=$( calc $n_threads/$n_cpu) ;

    cpu_usage=$( printf "%.0f" $cpu_usage )
    
    # echo "n_cpu" $n_cpu ", n_threads" $n_threads ", usage" $cpu_usage "%" ;
    echo $cpu_usage ;
}

cpu_usage=$(cpu_use)
# echo " cpu_usage" $cpu_usage "%"
 
if [ $cpu_usage -gt 90 ]; then
    echo " CPU_USAGE > 90.0 %, sleeping for a while ..." 
fi

while [ $cpu_usage -gt 90 ]; 
do
    sleep 5s ;
    cpu_usage=$(cpu_use)    
done
