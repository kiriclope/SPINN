#!/bin/bash

LC_NUMERIC=C
mem_usage=$(free -m | awk 'NR==2{print $3/$2*100}')
mem_usage=$(printf "%.0f" $mem_usage)

if [ $mem_usage -gt 85 ]; then
    echo " MEM_USAGE > 85.0%, sleeping for a while ..."
fi

while [ $mem_usage -gt 85 ]; do
    sleep 10s
    
    mem_usage=$(free -m | awk 'NR==2{print $3/$2*100}')
    mem_usage=$(printf "%.0f" $mem_usage)
done
