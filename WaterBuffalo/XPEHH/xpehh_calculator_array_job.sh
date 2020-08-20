#!/bin/bash
# Grid Engine options
#$ -N xpehh
#$ -cwd
#$ -M prasundutta87@gmail.com
#$ -m bea
#$ -pe sharedmem 4
#$ -l h_vmem=10G
#$ -t 1:24
#$ -l h_rt=24:00:0
#$ -R y
# Initialise the modules framework
. /etc/profile.d/modules.sh

chrom=`head -$SGE_TASK_ID chromosome_numbers.txt | tail -1`

breeds=('Jaffrabadi' 'Bhadawari' 'Jaffrabadi' 'Murrah' 'Jaffrabadi' 'Surti' 'Jaffrabadi' 'Banni' 'Jaffrabadi' 'Pandharpuri' 'Jaffrabadi' 'Med' 'Bhadawari' 'Murrah' 'Bhadawari' 'Surti' 'Bhadawari' 'Banni' 'Bhadawari' 'Pandharpuri' 'Bhadawari' 'Med' 'Murrah' 'Surti' 'Murrah' 'Banni' 'Murrah' 'Pandharpuri' 'Murrah' 'Med' 'Surti' 'Banni' 'Surti' 'Pandharpuri' 'Surti' 'Med' 'Banni' 'Pandharpuri' 'Banni' 'Med' 'Pandharpuri' 'Med')

for((j=0;j<=40;j=j+2))
do
	breed1=${breeds[j]}
	breed2=${breeds[j+1]}
	mkdir -p "$breed1"_vs_"$breed2"
	../hapbin-1.3.0/build/xpehhbin --hapA "$breed1"_all_chrom_hapbin_inputs/"$breed1"_hapbin_input_chr_"$chrom".impute.hap --hapB "$breed2"_all_chrom_hapbin_inputs/"$breed2"_hapbin_input_chr_"$chrom".impute.hap --map chr_"$chrom".map --out "$breed1"_vs_"$breed2"_chr_"$chrom"

mv "$breed1"_vs_"$breed2"_chr_"$chrom" "$breed1"_vs_"$breed2"
done





