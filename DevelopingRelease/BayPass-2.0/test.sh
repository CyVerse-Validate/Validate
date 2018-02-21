#!/bin/bash

#SBATCH -p development
#SBATch -t 00:30:00
#SBATCH -n 15
#SBATch -A iPlant-Collabs
#SBATCH -J test-baypass

inputGfile="../examples/lsa.geno"

tar -zxvf baypass_2.1.tar.gz
cd baypass_2.1_2/sources

make clean all FC=ifort

./i_baypass -npop 1 -gfile $inputGfile
