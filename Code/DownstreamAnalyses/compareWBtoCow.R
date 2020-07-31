library(tidyverse)
#library(data.table)
library(GenomicRanges)
#library(rtracklayer)

numPerms<-100

args = commandArgs(trailingOnly=TRUE)

#the parameters for peak calling
xpclr_high_Cow<-args[3]
xpclr_low_Cow<-args[4]
xpehh_high_Cow<-args[3]
xpehh_low_Cow<-args[4]

xpclr_high_WB<-args[1]
xpclr_low_WB<-args[2]
xpehh_high_WB<-args[1]
xpehh_low_WB<-args[2]

smoothed<-TRUE

#species<-"WaterBuffalo"

#setwd("C:/Users/jgdpr/Dropbox/Analysis/WaterBuffalo/")
setwd("/home/jprende3/Scratch/Andrea/")

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


liftOver<-function(param)
{
  lifted<-read_tsv(paste("peaks_",param,"_toHereford.bed", sep=""), col_names = c("Chr_COW", "Start_COW", "End_COW", "ID"))
  lifted<-lifted %>% separate(ID, into=c("Chr_WB", "Start_WB", "End_WB", "Contig_end"), convert = TRUE) %>%
    unite("ID", Chr_WB:End_WB, remove = FALSE)
  lifted$Width_WB<-lifted$End_WB-lifted$Start_WB
  lifted<-lifted %>% select(Chr_COW, End_COW, ID, Contig_end, Width_WB) %>% spread(Contig_end, End_COW)
  #sometimes both ends of peak could not be lifted over so remove these
  lifted<-lifted[complete.cases(lifted),]
  lifted$Width_Cow<-abs(lifted$END-lifted$START)
  lifted<-lifted %>%rowwise() %>% mutate(minC=min(END, START)) %>% mutate(maxC=max(END, START))
  #remove those regions where twice as big or small after lifting over
  lifted<-lifted[which(abs(log2(lifted$Width_WB/lifted$Width_Cow))<1),]
  lifted<-lifted %>% select(Chr=Chr_COW, start=minC, end=maxC, id=ID)
  return(makeGRangesFromDataFrame(lifted, keep.extra.columns=TRUE))
}

performPerm<-function(gr_WB, gr_Cow)
{
  #randomly select regions from the genome of the same size distribution of the XPCLR regions
  
  perm_list<-list()
  for(n in 1:length(gr_WB))
  {
    print(n)
    regSize<-end(gr_WB)[n]-start(gr_WB)[n]
    #for(i in 1:3)
    #{
    perm_list[[n]]<-generateRandomPos(numPerms,1:29, chroms_DF$Size, regSize)
    #}
  }
  
  #convert the permutations to grange objects and redo overlap with real regions to see if generally less or not
  permOverlaps<-list()
  for(m in 1:numPerms)
  {
    gr_perm<-makeGRangesFromDataFrame(as_tibble(t(sapply(perm_list,'[',c(m,m+numPerms,m+(2*numPerms)))), .name_repair = ~ c("chrom", "start", "end")))
    seqlevels(gr_perm)<-as.character(1:29)
    overlaps<-findOverlaps(gr_perm, gr_Cow)
    #overlaps
    permOverlaps[[paste("perm",m, sep="")]]<-c(m, length(unique(queryHits(overlaps))), length(unique(subjectHits(overlaps))))
  }
  
  permCounts<-as_tibble(t(bind_rows(permOverlaps)), .name_repair = ~ c("perm", "WB_count", "Cow_count"))
  overlaps<-findOverlaps(gr_WB, gr_Cow)
  permMetrics<-permCounts %>% summarise(WB_mean=mean(WB_count), WB_sd=sd(WB_count),Cow_mean=mean(Cow_count), Cow_sd=sd(Cow_count))
  permMetrics$numPerm<-numPerms
  permMetrics$WB_real<-length(unique(queryHits(overlaps)))
  permMetrics$Cow_real<-length(unique(subjectHits(overlaps)))
  permMetrics$WB_totalNum<-length(gr_WB)
  permMetrics$Cow_totalNum<-length(gr_Cow)
  permMetrics$WB_totalWidth<-sum(width(gr_WB))
  permMetrics$Cow_totalWidth<-sum(width(gr_Cow))
  permMetrics$WB_z<-(permMetrics$WB_real-permMetrics$WB_mean)/permMetrics$WB_sd
  permMetrics$Cow_z<-(permMetrics$Cow_real-permMetrics$Cow_mean)/permMetrics$Cow_sd
  permMetrics$WB_p<-2*pnorm(-abs(permMetrics$WB_z))
  permMetrics$Cow_p<-2*pnorm(-abs(permMetrics$Cow_z))
  return(permMetrics)
}

#load WB lifted to cattle
param_XPEHH<-paste(xpehh_high_WB, xpehh_low_WB, "WaterBuffalo_XPEHH", smoothed, sep="_")
wbOnCow_XPEHH_gr<-liftOver(param_XPEHH)
#restrict to peaks lifted to autosomes as some can be lifted to sex chroms or contigs
wbOnCow_XPEHH_gr<-wbOnCow_XPEHH_gr[which(seqnames(wbOnCow_XPEHH_gr) %in% 1:29)]
seqlevels(wbOnCow_XPEHH_gr)<-as.character(1:29)

param_XPCLR<-paste(xpclr_high_WB, xpclr_low_WB, "WaterBuffalo_XPCLR", smoothed, sep="_")
wbOnCow_XPCLR_gr<-liftOver(param_XPCLR)
wbOnCow_XPCLR_gr<-wbOnCow_XPCLR_gr[which(seqnames(wbOnCow_XPCLR_gr) %in% 1:29)]
seqlevels(wbOnCow_XPEHH_gr)<-as.character(1:29)


#load cattle peaks
load(paste("all_XPEHH_",xpehh_high_Cow, "_", xpehh_low_Cow, "_Cow_", smoothed, ".RData", sep=""))
gr_XPEHH<-reducePeaks(lowess_loci_dat_XPEHH)
#remove short XPEHH peaks to be more comparable to xpclr as some can be just one or two variants long
gr_XPEHH<-gr_XPEHH[which(width(gr_XPEHH)>5000),]
load(paste("all_XPCLR_",xpclr_high_Cow, "_", xpclr_low_Cow, "_Cow_", smoothed, ".RData", sep=""))
gr_XPCLR<-reducePeaks(lowess_loci_dat_XPCLR)

chroms_DF<-read_tsv("ARS-UCD1.2_chrom_sizes.txt", col_names = c("Chr", "Size"))
#restrict to autosomes
chroms_DF<-chroms_DF[1:29,]



metrics1<-performPerm(wbOnCow_XPEHH_gr, gr_XPEHH)
metrics1$xpehh_high_WB<-xpehh_high_WB
metrics1$xpehh_low_WB<-xpehh_low_WB
metrics1$xpclr_high_WB<-NA
metrics1$xpclr_low_WB<-NA
metrics1$xpehh_high_Cow<-xpehh_high_Cow
metrics1$xpehh_low_Cow<-xpehh_low_Cow
metrics1$xpclr_high_Cow<-NA
metrics1$xpclr_low_Cow<-NA

metrics2<-performPerm(wbOnCow_XPEHH_gr, gr_XPCLR)
metrics2$xpehh_high_WB<-xpehh_high_WB
metrics2$xpehh_low_WB<-xpehh_low_WB
metrics2$xpclr_high_WB<-NA
metrics2$xpclr_low_WB<-NA
metrics2$xpehh_high_Cow<-NA
metrics2$xpehh_low_Cow<-NA
metrics2$xpclr_high_Cow<-xpclr_high_Cow
metrics2$xpclr_low_Cow<-xpclr_low_Cow

metrics3<-performPerm(wbOnCow_XPCLR_gr, gr_XPEHH)
metrics3$xpehh_high_WB<-NA
metrics3$xpehh_low_WB<-NA
metrics3$xpclr_high_WB<-xpclr_high_WB
metrics3$xpclr_low_WB<-xpclr_low_WB
metrics3$xpehh_high_Cow<-xpehh_high_Cow
metrics3$xpehh_low_Cow<-xpehh_low_Cow
metrics3$xpclr_high_Cow<-NA
metrics3$xpclr_low_Cow<-NA

metrics4<-performPerm(wbOnCow_XPCLR_gr, gr_XPCLR)
metrics4$xpehh_high_WB<-NA
metrics4$xpehh_low_WB<-NA
metrics4$xpclr_high_WB<-xpclr_high_WB
metrics4$xpclr_low_WB<-xpclr_low_WB
metrics4$xpehh_high_Cow<-NA
metrics4$xpehh_low_Cow<-NA
metrics4$xpclr_high_Cow<-xpclr_high_Cow
metrics4$xpclr_low_Cow<-xpclr_low_Cow

metrics<-bind_rows(metrics1, metrics2, metrics3, metrics4)
write_tsv(metrics, paste("crossSpeciesPermMetrics_",xpehh_high_WB, "_", xpehh_low_WB, "_",xpclr_high_WB, "_", xpclr_low_WB, "_WBtoCow_",xpehh_high_Cow, "_", xpehh_low_Cow, "_",xpclr_high_Cow, "_", xpclr_low_Cow, "_Cow.txt", sep=""))


