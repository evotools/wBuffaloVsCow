# Plot the treemix results
args = commandArgs(TRUE)
inputfile = args[1]
popord = args[2]
outputname = args[3]
sourcedata = paste(args[4],'plotting_funcs.R', sep='/')
source(sourcedata)

# Plot treemix
pdf(paste(outputname, '_treemix.pdf', sep = ''))
plot_tree(inputfile, cex = 0.8, ybar = 2)
dev.off()
# Plot residuals
pdf(paste(outputname, '_residual.pdf', sep = ''))
plot_resid(inputfile, popord)
dev.off()

