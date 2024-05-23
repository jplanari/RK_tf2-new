#!/bin/bash

simfile=$1
rayleigh=$2
prandtl=$3
laplacian_discretization=$4
rkMethod=$5
rkFct=$6
nSims=$7

export GENERAL_MODS=$HOME/JOSEP/RK_tf2_mn5/mods
export SPECIFIC_MODS=$PWD/specific_mods
export CURR_DIR=$PWD

if [ "$#" -ne 7 ]; then
   echo "./launcher.sh <simfile> <rayleigh> <prandtl> <laplacian_discretization> <rkMethod> <rkFct> <nSims>"
   exit 0
fi

projectname="upc61"
bsc_scratch="/gpfs/scratch/${projectname}/JOSEP/DHC"
qosname="gp_resa"
jobname="PIT_DHC_${laplacian_discretization}_${pcgIters}_${rkMethod}_${rayleigh}"

dirname="${laplacian_discretization}_Ra${rayleigh}_Pr${prandtl}_${rkMethod}_${rkFct}_${nSims}"
foldername="${bsc_scratch}/${dirname}"
cfgname="${dirname}.cfg"

if ! test -d $foldername; then
   mkdir $foldername
fi

cp $simfile $foldername/$cfgname
cd $foldername

mkdir data timeav probes

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

echo "Solver_HPC2:Pressure_Solver:CG:Jacobi:Vol_Lap_CC:1000:1e-4" >> $cfgname

sbatch -A $projectname -q $qosname -J $jobname -t 20:00:00 $CURR_DIR/rs_launcher $cfgname 

