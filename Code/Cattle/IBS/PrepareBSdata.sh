#!/bin/bash
#$ -N prepibs             
#$ -cwd                  
#$ -l h_rt=23:59:59 
#$ -pe sharedmem 4
#$ -l h_vmem=16G
#$ -R y
#$ -e ./LOGS/$JOB_NAME.$JOB_ID.err
#$ -o ./LOGS/$JOB_NAME.$JOB_ID.out


# Initialise the environment modules
. /etc/profile.d/modules.sh
module load igmm/apps/plink/1.90b4
module load python/2.7.10


while getopts ":i:a:b:n" opt; do
  case $opt in
    i) infile=${OPTARG};;
    a) autosome=${OPTARG};;
    b) bootstr=${OPTARG};;
    n) nversion=${OPTARG};;
  esac
done

export PATH=/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/phylip/exe:$PATH

# Bootstrapped phylogeny tree.
plink --threads $NSLOTS --chr-set $autosome --maf 0.05 --bfile $infile --recode transpose --out ./ForBootstrappedDistances

python MakeBootstrapLists.py ./ForBootstrappedDistances $bootstr

