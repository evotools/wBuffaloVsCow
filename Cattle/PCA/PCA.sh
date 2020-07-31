#!/bin/bash
#$ -N pca
#$ -cwd                  
#$ -l h_rt=05:59:59 
#$ -pe sharedmem 4
#$ -l h_vmem=16G
#$ -R y

. /etc/profile.d/modules.sh
module load igmm/apps/plink/1.90b4
module load igmm/apps/vcftools
module load igmm/apps/tabix

vcftools --gzvcf $1 --max-alleles 2 --min-alleles 2 --recode --recode-INFO-all --stdout | bgzip -c > biallelic.vcf.gz && tabix -p vcf biallelic.vcf.gz


plink --cow --vcf biallelic.vcf.gz --vcf-min-qual 40 --update-ids updateIDs.txt --threads ${NSLOTS} --pca 410 --maf .05 --geno .2 --out PCA 
