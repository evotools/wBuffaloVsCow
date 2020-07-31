args = commandArgs(T)
library(tidyverse)

liftOver<-function(lifted, removed)
{
  # lifted<-read_tsv(infile, col_names = c("Chr_TGT", "Start_TGT", "End_TGT", "ID"))
  # lifted<-lifted %>% 
  #             separate(ID, into=c("XPCLR", "REGID"), convert = TRUE, sep = "#") %>% 
  #             separate(REGID, into=c("Chr_SRC", "Start_SRC", "End_SRC", "Contig_end"), convert = TRUE) %>%
  #             unite("ID", Chr_SRC:End_SRC, remove = FALSE)
  lifted = lifted[!lifted$ID %in% removed, ]
  lifted$Width_SRC<-lifted$End_SRC-lifted$Start_SRC
  lifted<-lifted %>% select(Chr_TGT, End_TGT, XPCLR, ID, Contig_end, Width_SRC) %>% spread(Contig_end, End_TGT)
  #sometimes both ends of peak could not be lifted over so remove these
  lifted<-lifted[complete.cases(lifted),]
  lifted$Width_TGT<-abs(lifted$END-lifted$START)
  lifted<-lifted %>%rowwise() %>% mutate(minC=min(END, START)) %>% mutate(maxC=max(END, START))
  return(lifted)
  
}

removeUnpaired <- function(file){
  starts = file2check %>% filter(Contig_end == "START")
  ends = file2check %>% filter(Contig_end == "END")
  starts = starts[starts$ID %in% ends$ID,] %>% select(ID, START_CHR = Chr_TGT)
  ends = ends[ends$ID %in% starts$ID,] %>% select(ID, END_CHR = Chr_TGT)
  joint = merge(starts, ends, by = "ID" )
  toRemove = joint[joint$START_CHR != joint$END_CHR, "ID"]
  return(toRemove)
}

fixOverlaps2 = function(lifted){
  for (n in c(2: nrow(lifted))){
    bpe = lifted[n-1,"maxC"]
    bpi = lifted[n,"minC"]
    if (bpe > bpi){
      lifted[n-1,"maxC"] = floor((bpe+bpi) / 2) - 1
      lifted[n,"minC"] = floor((bpe+bpi) / 2)
    }
  }
  return(lifted)
}

fixOverlaps <- function(lifted){
  overlapping <- which( (lifted[ 1:(nrow(lifted)-1),"Chr_TGT"] == lifted[ 2:(nrow(lifted)),"Chr_TGT"]) & ( lifted[ 1:(nrow(lifted)-1),"maxC"] > lifted[ 2:(nrow(lifted)),"minC"]) )
  midpoint = as.numeric(unlist(floor((lifted[overlapping, "maxC"] + lifted[overlapping+1, "minC"])/2)))
  lifted[overlapping, "maxC"] = midpoint - 1
  lifted[overlapping+1, "minC"] = midpoint 
  return(lifted)
}

for (mybed in args){
  file2check = read_tsv(mybed, col_names = c("Chr_TGT", "Start_TGT", "End_TGT", "ID"))%>% 
    separate(ID, into=c("XPCLR", "REGID"), convert = TRUE, sep = "#") %>% 
    separate(REGID, into=c("Chr_SRC", "Start_SRC", "End_SRC", "Contig_end"), convert = TRUE) %>%
    unite("ID", Chr_SRC:End_SRC, remove = FALSE)
  toRemove = removeUnpaired(file2check)
  lifted = liftOver(file2check, toRemove) 
  # Detect unusually big windows and remove them
  threshold = median(lifted$Width_TGT) + (5*IQR(lifted$Width_TGT))
  lifted = lifted %>% filter(Width_TGT <= threshold) %>% select(Chr_TGT, minC, maxC, XPCLR) %>% arrange(Chr_TGT, minC)
  # Fix eventual overlaps iteratively
  noverlaps = sum((lifted[ 1:(nrow(lifted)-1),"Chr_TGT"] == lifted[ 2:(nrow(lifted)),"Chr_TGT"]) & ( lifted[ 1:(nrow(lifted)-1),"maxC"] > lifted[ 2:(nrow(lifted)),"minC"]))
  while(noverlaps > 0){
    lifted = fixOverlaps(lifted) %>% arrange(Chr_TGT, minC)
    to_swap = lifted$minC > lifted$maxC
    lifted[to_swap, c("minC", "maxC")] = lifted[to_swap, c("maxC", "minC")]
    noverlaps = sum((lifted[ 1:(nrow(lifted)-1),"Chr_TGT"] == lifted[ 2:(nrow(lifted)),"Chr_TGT"]) & ( lifted[ 1:(nrow(lifted)-1),"maxC"] > lifted[ 2:(nrow(lifted)),"minC"]))
  }
  outname = gsub('.bed','.fixLift.bed',mybed)
  write.table(lifted, outname, col.names = F, row.names = F, sep = "\t", quote = F)
}


