library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)

args = commandArgs(trailingOnly=TRUE)

#species<-"Cow"
species<-args[7]


#the parameters for peak calling
xpclr_high<-args[1]
xpclr_low<-args[2]
xpehh_high<-args[3]
xpehh_low<-args[4]
xpclr_smooth<-args[5]
xpehh_smooth<-args[6]

#setwd("C:/Users/jgdpr/Dropbox/Analysis/WaterBuffalo/")
setwd("/home/jprende3/Scratch/Andrea/")
path_XPCLR<-paste("XPCLR/",species,"/", sep="")
path_XPEHH<-paste("XPEHH/",species,"/", sep="")

numChrom<-24
if(species=="Cow")
{
  numChrom<-29
}


#as peaks for different breeds can overlap use this code to find and merge overlapping peaks
#think this was also primarily plagiarised from Rebecca
reducePeaks<-function(lowess_loci_dat)
{
  g1 <- lowess_loci_dat %>% dplyr::rename(chr = Chr, start = MinPosition, end = MaxPosition) %>% dplyr::select(-id) %>% unique()
  
  gr <- makeGRangesFromDataFrame(g1, keep.extra.columns = TRUE)
  
  hits <- findOverlaps(gr, ignore.strand=TRUE, 
                       drop.self=FALSE, drop.redundant=FALSE)
  ovpairs <- Pairs(gr, gr, hits=hits)
  pint <- pintersect(ovpairs, ignore.strand=TRUE)
  mcols(pint) <- data.frame(first(ovpairs), second(ovpairs))
  rpint <- reduce(pint)
  
  compsInPeaks<-findOverlaps(rpint, gr)
  rpint2<-cbind.data.frame(as.data.frame(rpint[queryHits(compsInPeaks),]), as.data.frame(gr[subjectHits(compsInPeaks),])$Comp)
  
  rpint2<-aggregate(rpint2[6], rpint2[-6], 
                    FUN = function(X) paste(unique(X), collapse=", "))
  colnames(rpint2)[6]<-"Comp"
  rpint2<-makeGRangesFromDataFrame(rpint2, keep.extra.columns = TRUE)
  return(rpint2)
}

#code for randomly selecting regions of the genome
generateRandomPos <- function(n,chr,chr.sizes,width=1){
  random_chr <- sample(x=chr,size=n,prob=chr.sizes,replace=T)
  random_pos <- sapply(random_chr,function(chrTmp){sample(chr.sizes[chr==chrTmp]-width,1)}) 
  #res <- GRanges(random_chr,IRanges(random_pos,random_pos+width))
  res<-cbind(random_chr, random_pos, random_pos+width)
  return(res)
}

metrics<-list()

########################LOAD XPCLR AND XPEHH PEAKS####################################


load(paste("all_XPCLR_",xpclr_high, "_", xpclr_low, "_", species, "_", xpclr_smooth, ".RData", sep=""))

load(paste("all_XPEHH_",xpehh_high, "_", xpehh_low, "_", species, "_", xpehh_smooth, ".RData", sep=""))


#################OVERLAP XPCLR PEAKS WITH XPEHH PEAKS#################################


if(species=="Cow")
{
  chroms_DF<-read_tsv("ARS-UCD1.2_chrom_sizes.txt", col_names = c("Chr", "Size"))
  #restrict to autosomes
  chroms_DF<-chroms_DF[1:numChrom,]
} else {
  chroms_DF<-read_tsv("water_buffalo_re_arranged_chrom_ref_genome.fa.fai", col_names = c("Chr", "Size"))
  #restrict to autosomes
  chroms_DF<-chroms_DF[1:numChrom,]
}

gr_XPCLR<-reducePeaks(lowess_loci_dat_XPCLR)
gr_XPEHH<-reducePeaks(lowess_loci_dat_XPEHH)
#remove short XPEHH peaks to be more comparable to xpclr as some can be just one or two variants long
gr_XPEHH<-gr_XPEHH[which(width(gr_XPEHH)>5000),]

#randomly select regions from the genome of the same size distribution of the XPCLR regions
numPerms<-100
perm_list<-list()
for(n in 1:length(gr_XPCLR))
{
  print(n)
  regSize<-end(gr_XPCLR)[n]-start(gr_XPCLR)[n]
  #for(i in 1:3)
  #{
  perm_list[[n]]<-generateRandomPos(numPerms,1:numChrom, chroms_DF$Size, regSize)
  #}
  
}

#convert the permutations to grange objects and redo overlap with xpehh regions to see if generally less or not
permOverlaps<-list()
for(m in 1:numPerms)
{
  gr_XPCLR_perm<-makeGRangesFromDataFrame(as_tibble(t(sapply(perm_list,'[',c(m,m+numPerms,m+(2*numPerms)))), .name_repair = ~ c("chrom", "start", "end")))
  overlaps<-findOverlaps(gr_XPCLR_perm, gr_XPEHH)
  #overlaps
  permOverlaps[[paste("perm",m, sep="")]]<-c(m, length(unique(queryHits(overlaps))), length(unique(subjectHits(overlaps))))
}
permCounts<-as_tibble(t(bind_rows(permOverlaps)), .name_repair = ~ c("perm", "xpclr_count", "xpehh_count"))
overlaps<-findOverlaps(gr_XPCLR, gr_XPEHH)
permMetrics<-permCounts %>% summarise(xpclr_mean=mean(xpclr_count), xpclr_sd=sd(xpclr_count),xpehh_mean=mean(xpehh_count), xpehh_sd=sd(xpehh_count))
permMetrics$numPerm<-numPerms
permMetrics$xpclr_real<-length(unique(queryHits(overlaps)))
permMetrics$xpehh_real<-length(unique(subjectHits(overlaps)))
permMetrics$xpclr_totalNum<-length(gr_XPCLR)
permMetrics$xpehh_totalNum<-length(gr_XPEHH)
permMetrics$xpclr_totalWidth<-sum(width(gr_XPCLR))
permMetrics$xpehh_totalWidth<-sum(width(gr_XPEHH))
permMetrics$xpclr_z<-(permMetrics$xpclr_real-permMetrics$xpclr_mean)/permMetrics$xpclr_sd
permMetrics$xpehh_z<-(permMetrics$xpehh_real-permMetrics$xpehh_mean)/permMetrics$xpehh_sd
permMetrics$xpclr_p<-2*pnorm(-abs(permMetrics$xpclr_z))
permMetrics$xpehh_p<-2*pnorm(-abs(permMetrics$xpehh_z))
permMetrics$xpclr_peakHigh<-xpclr_high
permMetrics$xpclr_peakLow<-xpclr_low
permMetrics$xpehh_peakHigh<-xpehh_high
permMetrics$xpehh_peakLow<-xpehh_low
permMetrics$xpclr_smoothed<-xpclr_smooth
permMetrics$xpehh_smoothed<-xpehh_smooth
permMetrics$species<-species
metrics[[paste(xpclr_high, xpclr_low, xpehh_high, xpehh_low)]]<-permMetrics
#set real counts to have permutation id of 0
permCounts<-permCounts %>% add_row(perm=0, xpclr_count=length(unique(queryHits(overlaps))), xpehh_count=length(unique(subjectHits(overlaps))))
write_tsv(permCounts, paste("permCounts_",xpehh_high, "_", xpehh_low, "_",xpclr_high, "_", xpclr_low, "_", species, "_", xpehh_smooth, "_", xpclr_smooth, ".txt", sep=""))
#realCounts<-c("REAL", length(unique(queryHits(overlaps))), length(unique(subjectHits(overlaps))))
          

write_tsv(bind_rows(metrics), paste("permMetrics_",xpehh_high, "_", xpehh_low, "_",xpclr_high, "_", xpclr_low, "_", species, "_", xpehh_smooth, "_", xpclr_smooth, ".txt", sep=""))
