#!/bin/sh

source $HOME/JOSEP/tfa-whole/set_env.sh

simfile=$1
rayleigh=$2
prandtl=$3
laplacian_discretization=$4
rkMethod=$5
rkFct=$6
nSims=$7

if [ "$#" -ne 7 ]; then
   echo "./launcher.sh <simfile> <rayleigh> <prandtl> <laplacian_discretization> <rkMethod> <rkFct> <nSims>"
   exit 0
fi

projectname="upc61"
bsc_scratch="/gpfs/scratch/${projectname}/JOSEP/DHC"

dirname="${laplacian_discretization}_Ra${rayleigh}_Pr${prandtl}_${rkMethod}_${rkFct}_${nSims}"
foldername_probe="${bsc_scratch}/${dirname}/probes"
probe_fileName="${dirname}.tfpb"

cd foldername_probe

sbatch -A $projectname -q $qosname -J $jobname -t 20:00:00 -n1 tf-probe $probe_fileName 

mkdir datFiles && mv *dat datFiles

