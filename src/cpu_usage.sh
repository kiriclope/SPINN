#!/usr/bin/env bash

calc() { 
    awk "BEGIN { print $* * 100 }"; 
}

cpu_use() {
    n_cpu=$(nproc --all)
    n_threads=$(top -bn1 | grep "ni," | awk '{print $2}' | cut -d ',' -f1)

    if [[ "$n_threads" == "" ]]; then
        n_threads=0
    fi

    cpu_usage=$(calc $n_threads/$n_cpu)
    
    # Use LC_NUMERIC=C to ensure that printf recognizes the period as the decimal separator
    cpu_usage=$(LC_NUMERIC=C printf "%.0f" $cpu_usage)
    
    echo $cpu_usage
}

cpu_usage=$(cpu_use)

if [ $cpu_usage -gt 90 ]; then
    echo " CPU_USAGE > 90.0 %, sleeping for a while ..." 
fi

while [ $cpu_usage -gt 90 ]; do
    sleep 5s
    cpu_usage=$(cpu_use)    
done
