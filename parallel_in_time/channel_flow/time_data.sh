#!/bin/bash

file_input=$1
file_output=$2
nrhs=$3

if [ ! -f "$file_output" ]; then
   echo -e "#n\tt_ite\tt_spmv\tt_axpy\tt_axty\tt_dot\tpct_spmv" >> $file_output
fi

ite=$(grep 'iteration' $file_input | awk '{print $2, $7}')
spmv=$(grep 'kernel_spmv' $file_input | awk '{print $2, $7, $NF}' | sed -n '3p')
axpy=$(grep 'kernel_axpy' $file_input | awk '{print $2, $7}' | sed -n '2p')
axty=$(grep 'kernel_axty' $file_input | awk '{print $2, $7}' | sed -n '2p')
dot=$(grep 'kernel_sdot' $file_input | awk '{print $2, $7}')

t_ite=$(echo "$ite" | awk '{print $2/$1}')
t_spmv=$(echo "$spmv" | awk '{print $2/$1}')
t_axpy=$(echo "$axpy" | awk '{print $2/$1}')
t_axty=$(echo "$axty" | awk '{print $2/$1}')
t_dot=$(echo "$dot" | awk '{print $2/$1}')
pct_spmv=$(echo "$spmv" | awk '{print $3}')

echo -e "${nrhs}\t${t_ite}\t${t_spmv}\t${t_axpy}\t${t_axty}\t${t_dot}\t${pct_spmv}" >> $file_output
