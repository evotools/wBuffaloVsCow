#!/bin/bash
#$ -N ibsdist             
#$ -cwd                  
#$ -l h_rt=96:59:59 
#$ -pe sharedmem 4
#$ -l h_vmem=16G
#$ -R y
#$ -e ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -o ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.out

# Initialise the environment modules
. /etc/profile.d/modules.sh
module load igmm/apps/plink/1.90b4
module load python/2.7.10
module load R
module load igmm/compilers/gcc


while getopts ":i:a:n" opt; do
  case $opt in
    i) infile=${OPTARG};;
    a) autosome=${OPTARG};;
    n) nversion=${OPTARG};;
  esac
done

export PATH=/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/phylip/exe:$PATH

# Bootstrapped phylogeny tree.
python ./IBS_singleBS.py ./ForBootstrappedDistances $SGE_TASK_ID $autosome $NSLOTS
echo -e "n" | python ./Make_Phylip_input.py BS_$SGE_TASK_ID".mdist" $SGE_TASK_ID


