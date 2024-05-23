#!/bin/bash

simfile="dhc_template.cfg"
rayleigh="1e10"
prandtl="0.71"
laplacian_discretizations=("first" "cross" "cube")
pcgIters=$2
rkMethod="heunRK3"
rkFct="0.01"
nSims=("1" "2" "4" "8")
nRuns=$1

for laplacian_discretization in ${laplacian_discretizations[@]}; do

   processFile="${laplacian_discretization}_${pcgIters}_${rkMethod}.dat"

   for n in ${nSims[@]}; do
      for ((i = 1; i <= nRuns; i++)); do
         foldername="${laplacian_discretization}_Ra${rayleigh}_Pr${prandtl}_${pcgIters}_${rkMethod}_${rkFct}_${n}-${i}"
         ./time_data.sh "${foldername}/out/*out" ${processFile} $n
      done
   done

done
