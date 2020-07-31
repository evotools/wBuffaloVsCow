import sys
import subprocess as sbp
import random as rn

tfile = sys.argv[1]             # Input tped/tfam suffix.
nboot = int(sys.argv[2])       # number of bootstrap for analysis (100 recommended).
autosomes = sys.argv[3]
threads= sys.argv[4]


# Lists
inds = []
markers = []

# Store input tfam.
for line in open(tfile + '.tfam'):
	inds.append(line)

# Read input tped.
for line in open(tfile + '.tped'):
	markers.append(line)

# Number of SNP to choose.
nsnp = len(markers)

# Perform bootstrap.
mrklist = [int(i) for i in open("./LISTS/BS_{}.txt".format(nboot))]
otfam = open("BS_" + str(nboot) + '.tfam', 'w')
[otfam.write('%s' % i) for i in inds]
otfam.close()
otped = open("BS_" + str(nboot) + '.tped', 'w')
[otped.write("%s" % markers[mrk]) for mrk in mrklist]
otped.close()
# Reorder tped
p = sbp.Popen("Rscript arrange.R {}".format(nboot), shell = True)
p.wait()
cmd = "plink --chr-set {0} --threads {1} --allow-no-sex --nonfounders --tfile BS_{2} --distance 1-ibs flat-missing square --out BS_{2}"
p = sbp.Popen(cmd.format(autosomes, threads, nboot), shell = True)
#p = sbp.Popen("plink --chr-set " + autosomes + " --threads " + threads + " --allow-no-sex --nonfounders --tped BS_" + str(nboot) + ".tped --tfam BS_" + str(nboot) + ".tfam --distance 1-ibs flat-missing square --out BS_" + str(nboot), shell = True)
p.wait()
p = sbp.Popen("rm BS_{0}.tped BS_{0}.tfam".format(nboot), shell = True)
p.wait()
print "Bootstrapped IBS completed."
