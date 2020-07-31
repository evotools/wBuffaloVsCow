#!/bin/bash
#$ -N dataprep
#$ -cwd                  
#$ -l h_rt=05:59:59 
#$ -pe sharedmem 4
#$ -l h_vmem=16G
#$ -R y
#$ -P roslin_ctlgh

. /etc/profile.d/modules.sh
module load igmm/apps/plink/1.90b4
module load igmm/apps/vcftools

# Filter data
plink --cow --vcf $1 --make-bed --maf 0.05 --out unrelated_4treemix --threads 4
plink --cow --bfile unrelated_4treemix --keep IID2KeepTreemix.txt --make-bed --out unrelated_4treemix_min5ID
plink --cow --bfile unrelated_4treemix_min5ID --freq --family --out plink_freq
gzip plink_freq.frq.strat
python plink2treemix.py plink_freq.frq.strat.gz treemix_input.gz
awk '{print $1}' unrelated_4treemix_min5ID.fam > poporder

