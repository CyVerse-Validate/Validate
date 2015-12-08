#!/bin/bash
WRAPPERDIR=$( cd "$( dirname "$0" )" && pwd )
#SBATCH -p development
#SBATCH -t 00:30:00
#SBATCH -n 15
#SBATCH -A iPlant-Collabs 
#SBATCH -J test-bayesR
bfile="$WRAPPERDIR/test_data/simdata"
out="$WRAPPERDIR/test_data/output/out"
tar -xvf bin.tgz
bayesR -bfile $bfile -out $out
rm -rf bin

inputBED="test_data/simdata.bed"
inputBIM="test_data/simdata.bim"
inputFAM="test_data/simdata.fam"
output="bayesRTest.txt"

. ./wrapper.sh