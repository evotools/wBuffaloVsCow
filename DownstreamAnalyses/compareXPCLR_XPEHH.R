library(slider)
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)

args = commandArgs(trailingOnly=TRUE)

#species<-"Cow"
species<-args[6]
annotate<-TRUE
smoothed<-args[5]

#the parameters for peak calling
xpclr_high<-args[1]
xpclr_low<-args[2]
xpehh_high<-args[3]
xpehh_low<-args[4]


#setwd("C:/Users/jgdpr/Dropbox/Analysis/WaterBuffalo/")
setwd("/home/jprende3/Scratch/Andrea/")
path_XPCLR<-paste("XPCLR/",species,"/", sep="")
path_XPEHH<-paste("XPEHH/",species,"/", sep="")

numChrom<-24
if(species=="Cow")
{
  numChrom<-29
}

######################FUNCTIONS############################################
#function for peak calling based on parameters above
callPeaks<-function(xpclr_dat_lowess, comp, chr, max, min)
{
  #this peaks code was largely plagarised from rebecca
  #Peaks above max cutoff and give them a ID
  xpclr_dat_lowess1 <- xpclr_dat_lowess %>% 
    mutate(g35 = ifelse(Lowess_Abs_Std > max, 1,0),
           PeakID = rleid(g35)) 
  
  #Find maximum peak and its position
  max_peak_lowess <- xpclr_dat_lowess1 %>%
    filter(g35 == 1) %>%
    group_by(PeakID) %>%
    mutate(MaxPeak = max(Lowess_Abs_Std),
           PeakPosition = case_when(MaxPeak==Lowess_Abs_Std ~ Position)) %>%
    dplyr::select(PeakID, PeakPosition, MaxPeak) %>%
    na.omit()
  
  xpclr_dat_lowess2 <- xpclr_dat_lowess1 %>%
    left_join(max_peak_lowess)
  
  #xpehh_dat_lowess2 %>% filter(Lowess_Abs_Std_XPEHH > 3.5)
  
  # Find where starts and end position which is above min and includes position with XP_EHH above max
  xpclr_dat_lowess3 <- xpclr_dat_lowess2 %>% 
    mutate(g3 = ifelse(Lowess_Abs_Std > min, 1,0),
           PeakID_3 = rleid(g3)) %>%
    group_by(PeakID_3) %>% 
    mutate(MinPosition = min(Position),
           MaxPosition = max(Position),
           IncludesPeak = ifelse(MinPosition <= PeakPosition & MaxPosition >= PeakPosition, 1,0)) %>%
    ungroup() %>%
    filter(IncludesPeak == 1) %>%
    dplyr::select(MinPosition, MaxPosition, PeakPosition, MaxPeak) %>%
    unique()
  
  #There are some regions which include multiple points so mergeing those
  xpclr_dat_lowess4 <- xpclr_dat_lowess3 %>% 
    mutate(min_max_position = paste(MinPosition, MaxPosition, sep="_")) %>%
    group_by(min_max_position) %>%
    mutate(MaxPeak2 = max(MaxPeak)) %>%
    filter(MaxPeak == MaxPeak2) %>%
    ungroup %>%
    dplyr::select(-min_max_position, -MaxPeak2)
  
  #xpehh_dat_lowess4
  
  xpclr_dat_lowess5 <- xpclr_dat_lowess4 %>% 
    mutate(Chr = chr) %>%
    dplyr::select(Chr, MinPosition, MaxPosition) %>% 
    unique() %>%
    mutate(id = row_number(), 
           Comp=comp) 
  
  
  return(xpclr_dat_lowess5)
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



outputPeaks<-function(gr, overlaps, metric)
{
  if(metric == "XPEHH")
  {
    id<-paste(xpehh_high, xpehh_low, species, metric, smoothed, sep="_")
  } else if(metric == "XPCLR") {
    id<-paste(xpclr_high, xpclr_low, species, metric, smoothed, sep="_")
  }
  if(!file.exists(paste("peaks_", id, ".bed", sep="")))
  {
    
    #lift over start and end of peaks separately in case indels between species lead to splits in liftover
    gr_df<-as.data.frame(gr)
    bedStarts<-gr_df %>% select(seqnames,start=start, end=start)
    bedStarts$id<-paste(gr_df$seqnames, gr_df$start, gr_df$end, "START", sep="_")
    bedStarts$start<-bedStarts$start-1
    
    bedEnds<-gr_df %>% select(seqnames,start=end, end=end)
    bedEnds$id<-paste(gr_df$seqnames, gr_df$start, gr_df$end, "END", sep="_")
    bedEnds$start<-bedEnds$start-1
    write_tsv(bind_rows(bedStarts, bedEnds), paste("peaks_",id,".bed", sep=""), col_names=FALSE)
    
    #need liftover executable and chain file in same folder
    #command<-paste("./liftOver ", "peaks_",id,".bed liftover.chn peaks_",metrics,"_toHereford.bed peaks_",metrics,"_unliftedToHereford.bed",sep="")
    #print(command)
    system2("./liftOver", args=c(paste("peaks_",id,".bed", sep=""), "liftover.chn", paste("peaks_",id,"_toHereford.bed", sep=""), paste("peaks_",id,"_unliftedToHereford.bed", sep="")))
    
  }
  
  metrics<-paste(xpehh_high, xpehh_low, xpclr_high, xpclr_low, species, metric, smoothed, sep="_")
  write_tsv(left_join(as.data.frame(gr), overlaps), paste("Lowess_Loci_Locations_",metrics,".txt", sep=""))
  
  if(species == "Cow")
  {
    features<-import("Genomes/GCF_002263795.1_ARS-UCD1.2_genomic.gff.gz")
    features<-features[which(features$type == "gene"),]
    features<-renameSeqlevels(features,c(NC_037328.1="1", NC_037329.1="2", NC_037330.1="3", NC_037331.1="4", NC_037332.1="5", NC_037333.1="6", NC_037334.1="7", NC_037335.1="8", NC_037336.1="9", NC_037337.1="10", NC_037338.1="11", NC_037339.1="12", NC_037340.1="13", NC_037341.1="14", NC_037342.1="15", NC_037343.1="16", NC_037344.1="17", NC_037345.1="18", NC_037346.1="19", NC_037347.1="20", NC_037348.1="21", NC_037349.1="22", NC_037350.1="23", NC_037351.1="24", NC_037352.1="25", NC_037353.1="26", NC_037354.1="27", NC_037355.1="28", NC_037356.1="29", NC_037357.1="X", NC_006853.1="MT"))
    
  } else {
    features<-import("Genomes/GCF_003121395.1_ASM312139v1_genomic.gff.gz")
    features<-features[which(features$type == "gene"),]
    features<-renameSeqlevels(features,c(NC_037545.1="1", NC_037546.1="2", NC_037547.1="3", NC_037548.1="4", NC_037549.1="5", NC_037550.1="6", NC_037551.1="7", NC_037552.1="8", NC_037553.1="9", NC_037554.1="10", NC_037555.1="11", NC_037556.1="12", NC_037557.1="13", NC_037558.1="14", NC_037559.1="15", NC_037560.1="16", NC_037561.1="17", NC_037562.1="18", NC_037563.1="19", NC_037564.1="20", NC_037565.1="21", NC_037566.1="22", NC_037567.1="23", NC_037568.1="24", NC_037569.1="X", NC_006295.1="MT"))
  }
  
  genes<-findOverlaps(gr, features)
  regions<-gr[queryHits(genes),]
  geneIDs<-features[subjectHits(genes),]$Name
  genes<-as_tibble(cbind.data.frame(regions, geneIDs))
  genes<-full_join(overlaps,genes)
  genes<-genes %>% mutate(Overlap=replace_na(Overlap, "NO"))
  write_tsv(genes, paste(metrics,"_GenesList.txt", sep=""))
  
  #note if region didnt overlap any genes they wont be represented in this set
  genes<-aggregate(genes[8], genes[-8], 
                   FUN = function(X) paste(unique(X), collapse=", "))
  write_tsv(genes, paste(metrics,"_GenesByRegion.txt", sep=""))
  
  write_tsv(as.data.frame(features$Name), paste(species, "Background_Genes.txt", sep="_"))
  
}

########################XPCLR####################################

#load raw xpclr output files unless output already exists
#then run peak calling on lowess smoothed data
if(file.exists(paste("all_XPCLR_",xpclr_high, "_", xpclr_low, "_", species, "_", smoothed,".RData", sep="")))
{
  load(paste("all_XPCLR_",xpclr_high, "_", xpclr_low, "_", species, "_", smoothed,".RData", sep=""))
} else {
  files<-list.files(path=path_XPCLR,pattern = ".xpclr$", recursive = TRUE)
  
  if(species=="Cow")
  {
    #restrict to those breeds with most samples
    groups<-c("Holstein", "Boran", "NDama", "Kenana", "Ogaden", "Ankole", "Kazakh", "Lingnan", "TibetanYellow", "Baoule", "Muturu", "Brahman", "IranianBosTaurus")
    breedPairs<-as_tibble(files) %>% separate(value, into=c("Breed1", "Breed2", NA, NA, "Chrom", NA))
    files<-files[which((breedPairs$Breed1 %in% groups) & (breedPairs$Breed2 %in% groups))]
  }
  
  xpclrL<-list()
  
  lowess_loci_dat_XPCLR <- NULL
  for(file in files)
  {
    
    dat<-read_tsv(paste(path_XPCLR,file, sep=""))
    dat$position<-dat$stop-25000
    for(j in 1:numChrom)
    {
      dat2<-dat %>% filter(chrom == j)
      #dont think not smoothing XPCLR is going to work but put in for completeness
      if(smoothed==TRUE)
      {
        #xpclr_dat <- data.frame(lowess(dat$xpclr_norm ~ dat$position,f = 0.005))
        #xpclr_dat <- xpclr_dat %>%
        #  dplyr::rename( Position = x, Lowess_Abs_Std =y)
        
        dat2$rollMean<-slide_vec(abs(dat2$xpclr_norm), mean, .before = 3, .after=3)
        xpclr_dat<-select(dat2, position, rollMean)
        names(xpclr_dat)<-c("Position", "Lowess_Abs_Std")
        
      } else {
        xpclr_dat <- dat2 %>% select(position, xpclr_norm) %>% dplyr::rename( Position = position, Lowess_Abs_Std  = xpclr_norm)
      }
      
      thesePeaks<-callPeaks(xpclr_dat, file, dat2$chrom[1], xpclr_high, xpclr_low)
      lowess_loci_dat_XPCLR <- rbind(lowess_loci_dat_XPCLR, thesePeaks)
    }
    
    xpclrL[[file]]<-dat
  }
  
  lowess_loci_dat_XPCLR<-lowess_loci_dat_XPCLR %>% separate(Comp, into=c("Comp", NA), sep="/")
  #lowess_loci_dat_XPCLR<-aggregate(lowess_loci_dat_XPCLR[4], lowess_loci_dat_XPCLR[-4], 
  #          FUN = function(X) paste(unique(X), collapse=", "))
  
  xpclr<-bind_rows(xpclrL, .id = "file" )
  
  xpclr<-xpclr %>% separate(file, into=c("breed1", "breed2", NA, NA, NA))
  xpclr$comparison <- paste(xpclr$breed1, xpclr$breed2)
  
  save(lowess_loci_dat_XPCLR,file=paste("all_XPCLR_",xpclr_high, "_", xpclr_low, "_", species, "_", smoothed,".RData", sep=""))
}


################XPEHHH#########################################

#used !!!!!!!!!!!!!!.R to convert raw hapbin output to more manageable RData files
if(file.exists(paste("all_XPEHH_",xpehh_high, "_", xpehh_low, "_", species, "_", smoothed,".RData", sep="")))
{
  load(paste("all_XPEHH_",xpehh_high, "_", xpehh_low, "_", species, "_", smoothed,".RData", sep=""))
} else {
  
  lowess_loci_dat_XPEHH <- NULL
  
  for(i in 1:numChrom)
  {
    print(i)
    load(paste(path_XPEHH, "XPEHH_list.",i,".RData", sep=""))
    allComp$Pos <- as.numeric(allComp$Pos)
    for(breeds in unique(allComp$Comp))
    {
      #breeds<-"Bhadawari_Murrah"
      print(breeds)
      dat <- allComp %>% filter(Comp == breeds)
      
      #allComp<-allComp %>% separate(Comp,into=c("Pop1", "Pop2"), sep="_")
      
      if(smoothed == TRUE)
      {
        
        #this approach for smoothing first but seems less necessary for XPEHH
        #xpehh <-data.frame(lowess(abs(dat$`std XPEHH`) ~ dat$Pos,f = 0.000005))
        dat$rollMean<-slide_vec(abs(dat$`std XPEHH`), mean, .before = 500, .after=500)
        xpehh<-select(dat, Pos, rollMean)
        names(xpehh)<-c("Position", "Lowess_Abs_Std")
        #xpehh <- xpehh %>%
        #  dplyr::rename( Position = x, Lowess_Abs_Std =y)
        #max(xpehh$Lowess_Abs_Std)
        
        #ggplot(dat, aes(Pos, abs(`std XPEHH`)))+geom_point() + geom_line(aes(Pos, rollMean, colour="red"))
        #thesePeaks<-callPeaks(xpehh_dat_lowess, breeds, i, 3, 2)
      } else {
        #this for on unsmoothed data
        xpehh<-select(dat, Pos, `std XPEHH`)
        names(xpehh)<-c("Position", "Lowess_Abs_Std")
        xpehh$Lowess_Abs_Std<-abs(xpehh$Lowess_Abs_Std)
      }
      thesePeaks<-callPeaks(xpehh, breeds, i, xpehh_high, xpehh_low)
      lowess_loci_dat_XPEHH <- rbind(lowess_loci_dat_XPEHH, thesePeaks)
      
    }
  }
  save(lowess_loci_dat_XPEHH, file=paste("all_XPEHH_",xpehh_high, "_", xpehh_low, "_", species, "_", smoothed,".RData", sep=""))
}


#################OVERLAP XPCLR PEAKS WITH XPEHH PEAKS##############


if(annotate==TRUE)
{
  
  if(species=="Cow")
  {
    
  } else {
    chroms_DF<-read_tsv("water_buffalo_re_arranged_chrom_ref_genome.fa.fai", col_names = c("Chr", "Size"))
    #restrict to autosomes
    chroms_DF<-chroms_DF[1:numChrom,]
  }
  
  gr_XPCLR<-reducePeaks(lowess_loci_dat_XPCLR)
  gr_XPEHH<-reducePeaks(lowess_loci_dat_XPEHH)
  #remove short XPEHH peaks to be more comparable to xpclr as some can be just one or two variants long
  gr_XPEHH<-gr_XPEHH[which(width(gr_XPEHH)>5000),]
  overlaps<-findOverlaps(gr_XPCLR, gr_XPEHH)
  
  
  XPCLR_overlaps<-as.data.frame(gr_XPCLR[unique(queryHits(overlaps)),])
  XPCLR_overlaps$Overlap<-"YES"
  XPEHH_overlaps<-as.data.frame(gr_XPEHH[unique(subjectHits(overlaps)),])
  XPEHH_overlaps$Overlap<-"YES"
  
  #subsetByOverlaps(gr_XPEHH, gr_XPCLR)
  sum(width(subsetByOverlaps(gr_XPEHH, gr_XPCLR)))
  sum(width(subsetByOverlaps(gr_XPCLR, gr_XPEHH)))
  
  #param_XPCLR<-paste(xpehh_high, xpehh_low, xpclr_high, xpclr_low, species, "XPCLR", sep="_")
  #param_XPEHH<-paste(xpehh_high, xpehh_low, xpclr_high, xpclr_low, species, "XPEHH", sep="_")
  outputPeaks(gr_XPCLR, XPCLR_overlaps, "XPCLR")
  outputPeaks(gr_XPEHH, XPEHH_overlaps, "XPEHH")
  
}


