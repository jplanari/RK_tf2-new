#!/bin/bash

if [ "$#" -ne 2 ]; then
   echo "Error: ./nSIMS_RUN.sh <laplacian_discretization> <nRuns>"
   exit 0
fi

simfile="cf_template.cfg"
retau="320"
laplacian_discretization=$1
pcgIters="350"
rkMethod="heunRK3"
rkFct="0.8"
nSims=("1" "2" "4" "8")
nRuns=$2

if [ "$1" != "first" ] && [ "$1" != "cube" ] && [ "$1" != "cross" ]; then
   echo "Error: laplacian_discretization should be first, cube or cross. Please select correct input"
   exit 0
fi

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "SIMULATION PARAMETERS"
echo "      Frictional Reynolds number: ${retau}"
echo "      Discretization: ${laplacian_discretization}"
echo "      Maximum CG iterations: ${pcgIters}"
echo "      Runge-Kutta scheme: ${rkMethod}"
echo "      Safety factor: ${rkFct}"
echo "      Runs per nrhs: ${nRuns}"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

for n in ${nSims[@]}; do
   for ((i = 1; i <= nRuns; i++)); do
      echo "Launching case $i/$nRuns, nrhs=${n}"
      ./launcher.sh $simfile $retau $laplacian_discretization $pcgIters $rkMethod $rkFct $n $i
   done
done
