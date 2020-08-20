#!/bin/bash
# Grid Engine options
#$ -N hapbin_inputmaker_ne_358
#$ -cwd
#$ -M prasundutta87@gmail.com
#$ -m bea
#$ -pe sharedmem 1
#$ -l h_vmem=5G
#$ -t 1:7
#$ -l h_rt=2:00:0
#$ -R y
# Initialise the modules framework
. /etc/profile.d/modules.sh

module load igmm/apps/bcftools/1.6 igmm/apps/vcftools/0.1.13

sample=`head -$SGE_TASK_ID breeds.txt | tail -1`

bcftools view -S "$sample".txt beagle_output_phased_imputed_ne_358_new_sample_names_snp_id.vcf.gz -O z -o "$sample"_beagle_output_phased_imputed_ne_358.vcf.gz

for i in {1..24}
do
	vcftools --gzvcf "$sample"_beagle_output_phased_imputed_ne_358.vcf.gz --chr $i --IMPUTE --out "$sample"_hapbin_input_chr_$i
done

mkdir "$sample"_all_chrom_hapbin_inputs

mv "$sample"_hapbin_input* "$sample"_all_chrom_hapbin_inputs



