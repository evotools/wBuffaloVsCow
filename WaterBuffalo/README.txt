## Principal Component Analysis (PCA)

The PCAs were calculated using PLINK v1.90b4 64-bit and the following parameters:
  --allow-extra-chr
  --allow-no-sex
  --chr 1-24
  --chr-set 24
  --double-id
  --geno 0
  --maf 0.05
  --missing
  --nonfounders
  --out PCA_new_20
  --pca 30 var-wts header tabs
  --recode vcf
  --vcf ../Final_buffalo_biallelic_variants_81_samples_filtered.vcf.gz
  --vcf-min-gq 20

  The resulting .eigenvec file generated from PLINK was used in PCA_script.R to generate the PCA plot in Supplementary Figure 2.

## Variant Calling and Filtration
The tools and parameters used to perform this methodology has been presented in a tabular format in Variant_Calling_and_filtering_commands.docx

## TreeMix analysis
For Treemix analysis, at first, the filtered biallelic Single Nucleotide Variants (SNVs) were processed using PLINK v1.90b4 64-bit and following parameters:

  --allow-extra-chr
  --allow-no-sex
  --chr 1-24
  --chr-set 24
  --double-id
  --geno 0
  --maf 0.05
  --make-bed
  --nonfounders
  --remove ../../Admixture/resequenced_samples.txt
  --update-ids updatedIDs.txt
  --vcf ../../Final_buffalo_biallelic_variants_81_samples_filtered_with_snp_id.vcf.gz
  --vcf-min-gq 20

  The resulting PLINK binary file was again processed using PLINK v1.90b4 64-bit for calculating family-wise allele frequencies using these parameters:
  --bfile plink
  --chr 1-24
  --chr-set 24
  --family
  --freq
  --out for_treemix
  
A breed-stratified allele frequency report was obtained for 79 individuals using the above codes.

The breed-stratified allele frequency report was converted to an input file suitable for the TreeMix program using the python script ‘plink2treemix.py’ (https://bitbucket.org/nygcresearch/treemix/downloads/)

Finally 'treemix_graphmaker.R' was used to create Figure 1B.

## Compute XP-EHH
Chromosome-wise input files required for XP-EHH calculations for each breed was created using 'sample_wise_hapbin_inputfile_maker.sh'.

For calculating XP-EHH scores for all chromosomes and for all breed combinations (21 unique combinations for 7 breeds in this study), 'xpehh_calculator_array_job.sh' was run.

Both scripts can be submitted using the qsub command. For example: qsub sample_wise_hapbin_inputfile_maker.sh

The following breeds were present in the breeds.txt file in sample_wise_hapbin_inputfile_maker.sh script:

Banni
Bhadawari
Jaffrabadi
Med
Murrah
Pandharpuri
Surti

## Compute XP-CLR
Variants with global minor allele frequency < 1% were discarded using vcftools with the following command:
    
    vcftools --gzvcf input.vcf.gz --maf 0.01 --recode --recode-INFO-all --stdout | bgzip -c > filt.vcf.gz 
    
Individuals IDs were then separated into multiple sub-lists, named after the breed of interest, and saved into
one subfolder called LISTS.
We then created the pairwise comparisons of breeds using the script MakePairs.py as follow:
    
    python MakePairs.py ./LISTS/ 5 24 > toProcess.txt

Where LISTS is the subfolder which contains the lists of individuals, 5 is the minimum number of individuals 
per breed to consider for the analysis and 24 is the number of autosomes.
The output is a list of pairwise comparison, repeated 24 times (one per chromosome).
We then submit the analysis to the scheduler using the following command:
    
    qsub -t 1-`wc -l toProcess.txt | cut -f 1` -tc 150 01-XPCLR.sh  -f filt.vcf.gz -l toProcess.txt -w 50000 -s 600
    
where 50000 is the windows size and 600 is the max number of markers per windows.
Finally, we gathered the results using the following command:
    
    ./02-Combine.sh 24
    
Where 24 is the number of chromosomes processed. The final outputs are collected in XPCLR_RESULTS in
pairs of breeds (XPCLR_RESULTS/Breed1_Breed2).
