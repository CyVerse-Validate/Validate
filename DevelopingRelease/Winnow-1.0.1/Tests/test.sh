#!/bin/bash
WRAPPERDIR=$( cd "$( dirname "$0" )" && pwd )
#SBATCH -p development
#SBATCH -t 00:30:00
#SBATCH -n 15
#SBATCH -A iPlant-Collabs
#SBATCH -J test-winnow
knowntruth="$WRAPPERDIR/kt.ote"
Folder="$WRAPPERDIR/OutputPLINK/"
SNP="SNP"
Score="P"
beta="BETA"
Pvaladjust="fdr_bh"
seper="whitespace"
kttype="OTE"
module load python

#Now to execute the main program
python ../bin/winnow.py --Folder $Folder --Class $knowntruth --Snp $SNP --Score $Score --beta $beta --kttype $kttype --pvaladjust $pvaladjust

