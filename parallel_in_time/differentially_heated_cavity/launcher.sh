#!/bin/bash

simfile=$1
rayleigh=$2
prandtl=$3
laplacian_discretization=$4
pcgIters=$5
rkMethod=$6
rkFct=$7
nSims=$8
i=$9

if [ "$#" -ne 9 ]; then
   echo "./launcher.sh <simfile> <rayleigh> <prandtl> <laplacian_discretization> <pcgIters> <rkMethod> <rkFct> <nSims> <counter>"
   exit 0
fi

foldername="${laplacian_discretization}_Ra${rayleigh}_Pr${prandtl}_${pcgIters}_${rkMethod}_${rkFct}_${nSims}-${i}"
cfgname="${foldername}.cfg"

if ! test -d $foldername; then
   mkdir $foldername
fi

cp $simfile $foldername/$cfgname
cd $foldername

echo "Config:NumSims:${nSims}" >> $cfgname
echo "Param:int:nSims:${nSims}" >> $cfgname

if [ "$laplacian_discretization" == "first" ]; then
   echo "MatrixDefault:Vol_Lap_CC" >> $cfgname
else
   echo "Pattern:CC:p_CC_${laplacian_discretization}:smesh-high-order.so" >> $cfgname
   echo "Kernel:k_lap:k_Vol_Lap_CC_${laplacian_discretization}:smesh-high-order.so" >> $cfgname
   echo "Matrix:Vol_Lap_CC:CC:k_lap" >> $cfgname
fi

echo "Param:double:Ra:${rayleigh}" >> $cfgname
echo "Param:double:Pr:${prandtl}" >> $cfgname
echo "Param:double:RKfct:${rkFct}" >> $cfgname
echo "Param:str:RKmethod:${rkMethod}" >> $cfgname

echo "Solver_HPC2:Pressure_Solver:CG:Jacobi:Vol_Lap_CC:${pcgIters}:1e-12" >> $cfgname

projectname="upc61"
qosname="gp_resa"
jobname="PIT_DHC_${laplacian_discretization}_${pcgIters}_${rkMethod}_${rayleigh}"

sbatch -A $projectname -q $qosname -J $jobname -t 20:00:00 ../rs_launcher $cfgname 

