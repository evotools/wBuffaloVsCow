# XP-CLR to BOMA
## Introduction
This readme explains how to convert the XPCLR regions to BOMA.

## Dependencies
The workflow relies on the following dependencies:
 1. [R](https://www.r-project.org/)
 2. [LiftOver](https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver)
 3. [bedGraphToBigWig](https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig)
 4. [bedSort](https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedSort)

## Procedure summary - Cattle
First, convert the single XPCLR values to bed, finding the central portion of each, using ```convertToBed``` 
that removes 10Kb from each flank. 
Then, the start and end positions for each region are saved as a separate
BED entry using ```awk``` (vertical bed), lifted using ```liftOver``` tool from UCSC and sorted with ```bedtools sort```.
Finally, the vertical lifted bed is combined together using ```liftXPCLR.R``` R script that:
 1. Remove the regions with the delimiting region on different chromosomes
 2. Reconstitutes the regions keeping initial and end positions lifted to the target
 3. Drop regions with initial and end positions too far apart (median + 5*IQR)
 4. Fix eventual overlaps in the lifted regions
The resulting bed files are then converted to bigWig using bedGraphToBigWig, dropping the regions without values (NA)
and providing the chromosomes' sizes.
The procedure can be performed through the following code:
```
for i in `ls | grep _cattle_XPCLR_03062020.xpclr`; do 
	bname=`basename -s ".xpclr" ${i}`
	cat $i | ../convertToBed - 10000 > ${bname}.raw.bed
	awk 'BEGIN{n=1}; {print $1,$2-1,$2,$4"#"$1"_"$2"_"$3"_START";print $1,$3-1,$3,$4"#"$1"_"$2"_"$3"_END"}' ${bname}.raw.bed | \
		liftOver stdin liftover_HerefordToBuffalo.chn stdout ${bname}.unmapped.vert.bed \
		| bedtools sort -i - > ${bname}.lifted.vert.bed;
	Rscript liftXPCLR.R ${bname}.lifted.vert.bed
	cat ${bname}.lifted.vert.bed | awk '$4!="NA"{print}' | bedSort stdin ${bname}.bed
	bedGraphToBigWig ${bname}.bed WaterBuffaloChromosomes.sizes ${bname}_buffalo_XPCLR_03062020.bw
	echo "Done $i"
	echo 
done
```

## Procedure summary - Water Buffalo
The procedure for the Water Buffalo is more streamlined, since it doesn't require the liftover 
and region fixing steps.
The code is as follow:
```
for i in `ls | grep _buffalo_XPCLR_03062020.xpclr`; do
        bname=`basename -s ".xpclr" ${i}`
        cat $i | ../convertToBed - 10000 | bedSort - ${bname}.bed
	bedGraphToBigWig ${bname}.bed WaterBuffaloChromosomes.sizes ${bname}_buffalo_XPCLR_03062020.bw
done
```
