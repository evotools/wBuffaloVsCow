import sys
import subprocess as sbp
import random as rn
import os

tfile = sys.argv[1]             # Input tped/tfam suffix.
nboots = int(sys.argv[2])       # number of bootstrap for analysis (100 recommended).


# Lists
inds = []
markers = []

# Read input tped.
nsnp = len([line for line in open(tfile + '.tped')])

# Make list directory
if not os.path.isdir("./LISTS"): os.mkdir("./LISTS")

# Perform bootstrap.
for boot in xrange(0, nboots):
	boot += 1
	otped = open("./LISTS/BS_" + str(boot) + '.txt', 'w')
	[otped.write( "{}\n".format(rn.randint(0, nsnp-1)) ) for i in xrange(0, nsnp)]
	otped.close()
print "Bootstrapped lists ready."
