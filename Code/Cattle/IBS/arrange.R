library(tidyverse)
args = commandArgs(T)
nboot = args[1]

tped = read_delim(paste("BS_", nboot, ".tped", sep = ''), delim = ' ', col_names = F) %>% arrange(X1, X4) 
sortped = tped[order(tped$X1, tped$X4, decreasing = F),]
write.table(sortped, paste("BS_", nboot, ".tped", sep = ''), sep = ' ', col.names = F, row.names = F, quote = F)


