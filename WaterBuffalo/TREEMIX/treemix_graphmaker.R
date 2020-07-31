##This script will plot TreeMix phylogeny graphs for a variable number 
##of migration events denoted by 'i'. Currently, only the first migration 
##event has been plotted. The 'for' loop can be edited for plotting
##graphs for n number of migration events

setwd("~/81_samples_genomic_variation/Treemix/treemix_new/")

source("plotting_funcs.R")

pdf('treemix_output_migration_model_1.pdf', width=19, height=11)
par(mfrow=c(3,2))
for(i in 1:1){
  plot_tree(paste0("treemix_output_migration_model_",i,"_LD"),cex = 0.7)
}
dev.off()
