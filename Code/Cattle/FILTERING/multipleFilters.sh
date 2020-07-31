#!/bin/bash
#$ -N filt
#$ -l h_vmem=8G
#$ -l h_rt=23:59:59
#$ -pe sharedmem 1
#$ -R y
#$ -cwd
#$ -e LOGS/$JOB_NAME.$TASK_ID.err
#$ -o LOGS/$JOB_NAME.$TASK_ID.out

# Load modules
. /etc/profile.d/modules.sh
module load igmm/apps/vcftools

# Get thresholds
CRV=$( sed "$SGE_TASK_ID q;d" $2 | awk '{print $1}' )
GQV=$( sed "$SGE_TASK_ID q;d" $2 | awk '{print $2}' )

if [ ! -e CR${CRV} ]; then mkdir CR${CRV}; fi 
if [ ! -e CR${CRV}/GQ${GQV} ]; then mkdir CR${CRV}/GQ${GQV}; fi

if [ ${file: -4} == ".txt" ]; then
    vcf-concat -f $1 | vcftools --vcf - --min-alleles 2 --max-alleles 2 --minGQ $GQV --max-missing $CRV --out ./CR${CRV}/GQ${GQV}/filtering_GQ${GQV}_CR${CRV}
elif [ ${file: -7} == ".vcf.gz" ]; then
    vcftools --gzvcf $1 --min-alleles 2 --max-alleles 2 --minGQ $GQV --max-missing $CRV --out ./CR${CRV}/GQ${GQV}/filtering_GQ${GQV}_CR${CRV}
else
    vcftools --vcf $1 --min-alleles 2 --max-alleles 2 --minGQ $GQV --max-missing $CRV --out ./CR${CRV}/GQ${GQV}/filtering_GQ${GQV}_CR${CRV}
fi
echo "Done GQ = $GQV, missingness = $CRV"
