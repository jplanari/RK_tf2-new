#!/bin/bash

simfile=$1
rayleigh=$2
prandtl=$3
nprocs=$4
rkMethod=$5
rkFct=$6
nSims=$7

if [ "$#" -ne 7 ]; then
   echo "./launcher.sh <simfile> <rayleigh> <prandtl> <nprocs> <rkMethod> <rkFct> <nSims>"
   exit 0
fi

cfgname=$simfile

echo "Config:NumSims:${nSims}" >> $cfgname
echo "Param:int:nSims:${nSims}" >> $cfgname

echo "Param:double:Ra:${rayleigh}" >> $cfgname
echo "Param:double:Pr:${prandtl}" >> $cfgname
echo "Param:double:RKfct:${rkFct}" >> $cfgname
echo "Param:str:RKmethod:${rkMethod}" >> $cfgname

mpirun -n $nprocs tf-sim $simfile

