#!/bin/bash

simfile=$1
reynolds=$2
rkMethod=$3
rkFct=$4
nx=$5
ny=$6
nz=$7
tpt=$8
cName=$9

if [ "$#" -ne 9 ]; then
   echo "./launcher.sh <simfile> <reynolds> <rkMethod> <rkFct> <nx> <ny> <nz> <threads_per_task> <cName>"
   exit 0
fi

foldername="${cName}_Re${reynolds}_${rkMethod}_${rkFct}_${nx}${ny}${nz}"
cfgname="${foldername}.cfg"

if ! test -d $foldername; then
   mkdir $foldername
fi

nt=$((nx * ny * nz))

cp $simfile $foldername/$cfgname
cd $foldername

echo "SMesh:Part:${nx}:${ny}:${nz}" >> $cfgname

echo "Param:double:Re_tau:${reynolds}" >> $cfgname
echo "Param:double:RKfct:${rkFct}" >> $cfgname
echo "Param:str:RKmethod:${rkMethod}" >> $cfgname

echo "HPC2Config:-profile:-overlap:-thr=${tpt}" >> $cfgname

projectname="upc61"
qosname="gp_resa" 
jobname="CF_${rkMethod}_Re${reynolds}"

sbatch -A $projectname -q $qosname -J $jobname -t 20:00:00 --ntasks=$nt ../rs_launcher $cfgname 

