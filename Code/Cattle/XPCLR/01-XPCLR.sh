#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N xpclr
#$ -cwd
#$ -l h_rt=23:00:00
#$ -pe sharedmem 1
#$ -R y
#$ -l h_vmem=16.0G
#$ -P roslin_ctlgh
#$ -o ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.err

. /etc/profile.d/modules.sh
module load roslin/bcftools
module load igmm/apps/vcftools
module load igmm/apps/tabix
module load anaconda
source activate xpclr

while getopts ":f:l:w:s:n" opt; do
  case $opt in
    f) filename=${OPTARG};;
    l) list=${OPTARG};;
    w) winsize=${OPTARG};;
    s) maxsnp=${OPTARG};;
    n) nversion=${OPTARG};;
  esac
done

# Get input VCF
tgt=$( sed "$SGE_TASK_ID q;d" ${list} | awk '{print $1}' )
p1=$( sed "$SGE_TASK_ID q;d" ${list} | awk '{print $2}' )
p2=$( sed "$SGE_TASK_ID q;d" ${list} | awk '{print $3}' )
echo "Chromosome: $tgt"
echo "Population 1: $p1"
echo "Population 2: $p2"
echo "Window size: $winsize"
echo "Maximum SNPs per window: $maxsnp"

if [ ! -e OUTPUTS ]; then mkdir OUTPUTS; fi
if [ ! -e OUTPUTS/${p1}_${p2} ]; then mkdir OUTPUTS/${p1}_${p2}; fi

echo "xpclr -F vcf -I $filename -Sa LISTS/${p1}.txt -Sb LISTS/${p2}.txt -O ./OUTPUTS/${p1}_${p2}/${p1}_${p2}.$tgt.xpclr -C $tgt --size $winsize --maxsnps $maxsnp"
xpclr -F vcf -I $filename -Sa LISTS/${p1}.txt -Sb LISTS/${p2}.txt -O ./OUTPUTS/${p1}_${p2}/${p1}_${p2}.$tgt.xpclr -C $tgt --size $winsize --maxsnps $maxsnp

