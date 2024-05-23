#!/bin/bash

file_input=$1
file_output=$2
nrhs=$3

if [ ! -f "$file_output" ]; then
   echo -e "#n\tt_ite\tt_spmv\tt_axpy\tt_axty\tt_dot\tT_poisson\tT_spmv\tpct_poisson" >> $file_output
fi

ite=$(grep 'iteration' $file_input | awk '{print $2, $5}')
spmv=$(grep 'kernel_spmv' $file_input | awk '{print $2, $5}' | sed -n '5p')
poisson=$(grep 'poisson' $file_input | awk '{print $5, $NF}')
axpy=$(grep 'kernel_axpy' $file_input | awk '{print $2, $5}' | sed -n '2p')
axty=$(grep 'kernel_axty' $file_input | awk '{print $2, $5}' | sed -n '2p')
dot=$(grep 'kernel_sdot' $file_input | awk '{print $2, $5}')


t_ite=$(echo "$ite" | awk '{print $2/$1}')
T_poisson=$(echo "$poisson" | awk '{print $1}')
t_spmv=$(echo "$spmv" | awk '{print $2/$1}')
T_spmv=$(echo "$spmv" | awk '{print $2}')
t_axpy=$(echo "$axpy" | awk '{print $2/$1}')
t_axty=$(echo "$axty" | awk '{print $2/$1}')
t_dot=$(echo "$dot" | awk '{print $2/$1}')
pct_poisson=$(echo "$poisson" | awk '{print $2}')

echo -e "${nrhs}\t${t_ite}\t${t_spmv}\t${t_axpy}\t${t_axty}\t${t_dot}\t${T_poisson}\t${T_spmv}\t${pct_poisson}" >> $file_output
