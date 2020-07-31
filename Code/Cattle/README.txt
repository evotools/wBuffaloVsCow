# Cattle analyses

## Post-VQSR analyses
After completing the variant quality score recalibration of the variants and performing hard-filtering of 
InDels in the dataset, the individuals with high missingness or related (--relatedness2 > 0.0625) were
excluded from subsequent analyses.

## Principal Component Analysis
Principal component analysis (PCA) has been performed using plink v1.90b4 calculating all components.
We filtered out genotypes with minor allele frequency < 5% and call rate < 80%, non-biallelic and with QUAL<40: 

    plink --cow --geno .2 --maf .05 --out PCA --pca 410 --threads 4 --update-ids updateIDs.txt --vcf biallelic.vcf.gz --vcf-min-qual 40

This analysis can be run on a SGE cluster environment using the script in the PCA folder, using the submission
command:

    qsub ./PCA/PCA.sh myvcf.vcf.gz

## Genotype refinement
To define the filtering criteria for genotype quality (GQ) and call rate (CR), we computed all combinations of 
a series of thresholds for both parameters, creating a curve chart that allowed us to identify the most
balanced combination of parameters. The software need to have access to vcftools (>=v0.1.13)
To do so, we computed the number of markers retained applying each combination of 
filtering using vcftools:

    vcftools --gzvcf $1 --min-alleles 2 --max-alleles 2 --minGQ $GQV --max-missing $CRV --out ./CR${CRV}/GQ${GQV}/filtering_GQ${GQV}_CR${CRV}

Where CRV is the CR and GQV is the GQ value.
In order to speed up all the analyses, we parallelized the whole process in a SGE environment with the script in the script
./FILTERIG/multipleFilters.sh. The thresholds' combination can be defined as pairs of call rates and GQ (see ./FILTERING/thresholds.txt)
and the script submitted as follow:

    qsub -t `wc -l ./FILTERING/thresholds.txt | awk '{print $1}'` ./FILTERIG/multipleFilters.sh [myvcf.vcf.gz/myvcf.vcf/mylistofvcf.txt] ./FILTERING/thresholds.txt

The resulting logs can be gathered using the following commands:

    while read p; do
        CRV=`echo $p | awk '{print $1}'`
        GQV=`echo $p | awk '{print $2}'`
        SNPN=`grep Sites CR${CRV}/GQ${GQV}/filtering_GQ${GQV}_CR${CRV}.log` | grep kept | cut -f 4 -d ' ' `
        echo $CRV $GQV $SNPN
    done < thresholds.txt > results.txt

The file results.txt will include the number of markers retained by each filtering.


## Treemix
Treemix analysis has been performed using the filtered dataset reduced to populations with at least 5 individuals.
To run this scripts, the system need to have availability of plink (v1.90), R, treemix, threepop and fourpop.
The preparation of the dataset can be accessed in the file TREEMIX/prepareData.sh, whereas the treemix instruction are 
shown in TREEMIX/TreemixEddie.sh.
The script can be run in a SGE environment through the commands:

    qsub TREEMIX/prepareData.sh myvcf.vcf.gz
    qsub -hold_jid dataprep TREEMIX/TreemixEddie.sh -i unrelated_4treemix_min5ID -n 0 -r ./TREEMIX -l n -s 2700000000
    qsub -hold_jid dataprep -t 1-10 TREEMIX/TreemixEddie.sh -i unrelated_4treemix_min5ID -r ./TREEMIX -l n -s 2700000000

The script will take care of generating all the outputs, including the graphs. It is possible to consider linkage disequilibrium
within the analysis but changin -l n to -l y. The script will take care of computing the density based on the map of markers 
provided and the length of the genome.


## Identity by State analysis.
We coputed a phylogeny tree for the different breeds using identity by state (IBS)-based distances (1-IBS) with 100 bootstraps.
The script to perform the analysis are in folder ./IBS, and need to have installed plink (v1.90), python 2.7, R, and phylip.
The analysis has been carried out using the same input of Treemix to leverage the populations with low sample size.
Data were first pre-processed using plink:

    plink --bfile unrelated_4treemix_min5ID --chr-set 29 --maf 0.05 --out ./ForBootstrappedDistances --recode transpose --threads 4

Then, we create the bootstrapped variants lists using the script in ./IBS:

    python MakeBootstrapLists.py ./ForBootstrappedDistances $bootstr

Where 'bootstr' is the number of bootstrap (i.e. 100). The files generated are then the inputs for the next script, which can be
submitted to an SGE cluster as follow: 

    qsub -t 1-100 IBS.sh -i ForBootstrappedDistances -a AUTOSOMES

This script will simply call the script IBS_singelBS.py, which extracts the variants listed for each iterations, including repetitions.
The same script then compute the distances through plink running the command:
    
    plink --chr-set AUTOSOMES --threads THREADS --allow-no-sex --nonfounders --tfile BS_N --distance 1-ibs flat-missing square --out BS_N

Where 'N' is the number of bootstrap to be executed.
Each bootstrap iteration can be run manually, one by one, using the command:

    python IBS_singleBS.py ForBootstrappedDistances N AUTOSOMES THREADS

The submitted script then convert the output from plink to phylip using the python script Make_Phylip_input.py as follow:

    echo -e "n" | python ./Make_Phylip_input.py BS_$SGE_TASK_ID".mdist" N

The initial echo -e "n" specify to the script that the the tree is unrooted, so it will not return the row numbers of the root breeds.
Finally, it is possible to create the newick tree by using the script PHYLIP.sh

    qsub PHYLIP.sh 

This will automatically combine all outputs from the bootstrap and process it through neighbors and consense tools.
It is also possible to specify the outgroup number (as row number of one individual belonging to the outgroup) using
-o NOG, where NOG is the outgroup row number. The resulting tree can be visualised with the most common tree visualization tools
accepting newick as input format. However, to reproduce the figure in the paper, it is needed to run the steps in the next section
of this readme.

## Create tree
To create the tree, it is necessary to have installed:

 1. FigTree[https://github.com/rambaut/figtree]
 2. forester[https://sites.google.com/site/cmzmasek/home/software/forester/phyloxml-converter], 
 3. Graphlan[https://bitbucket.org/nsegata/graphlan/wiki/Home]
 4. Python 2.7 with the packages colormap, xml and biopython.

The tree is processed using FigTree as follow:

   a. Generate a tree with the branches coloured, but in the same order and length
   b. Save it as a Nexus file (file > export tree > Nexus > flag all options)
   c. Generate a transformed tree (left menu > tree > flag transform branches)
   d. Save it as nexus tree as well
 
Then, open the two trees in forester.jar and export them as PhyloXML 
   (file > save tree as... > select PhyloXML)

The two phyloxml files are then joined using JoinPhyloXMLannotations.py:

    python JoinPhyloXMLannotations.py Tree1.xml Tree2.xml

The resulting phyloxml is then fixed for marker scale and labels size using the script FixGraphlanXml.py:

    python FixGraphlanXml.py final.xml 10 1 > toplot.xml

Finally, run graphlan as follow:

    graphlan.py toplot.xml mytree.png --dpi 300 --size 15

This will generate the resulting phylogenetic tree. 

## Compute XP-CLR
Variants with global minor allele frequency < 1% were discarded using vcftools with the following command:
    
    vcftools --gzvcf input.vcf.gz --maf 0.01 --recode --recode-INFO-all --stdout | bgzip -c > filt.vcf.gz 
    
Individuals IDs were then separated into multiple sub-lists, named after the breed of interest, and saved into
one subfolder called LISTS.
We then created the pairwise comparisons of breeds using the script MakePairs.py as follow:
    
    python MakePairs.py ./LISTS/ 5 29 > toProcess.txt

Where LISTS is the subfolder which contains the lists of individuals, 5 is the minimum number of individuals 
per breed to consider for the analysis and 29 is the number of autosomes.
The output is a list of pairwise comparison, repeated 29 times (one per chromosome).
We then submit the analysis to the scheduler using the following command:
    
    qsub -t 1-`wc -l toProcess.txt | cut -f 1` -tc 150 01-XPCLR.sh  -f filt.vcf.gz -l toProcess.txt -w 50000 -s 600
    
where 50000 is the windows size and 600 is the max number of markers per windows.
Finally, we gathered the results using the following command:
    
    ./02-Combine.sh 29
    
Where 29 is the number of chromosomes processed. The final outputs are collected in XPCLR_RESULTS in
pairs of breeds (XPCLR_RESULTS/Breed1_Breed2).
 
