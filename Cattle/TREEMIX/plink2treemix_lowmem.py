import sys, os, gzip, numpy

if len(sys.argv) < 3 or len(sys.argv) > 4:
	print "plink2treemix.py [gzipped input file] pops [gzipped output file]"
	print "ERROR: wrong command line"
	exit(1)
infile = gzip.open(sys.argv[1])
inpop = open(sys.argv[2])
try: autosome = sys.argv[4]
except: autosome = None
if autosome is None:
    outfile = gzip.open(sys.argv[3], "w")
else:
    outfile = gzip.open(sys.argv[3].replace(".gz", "."+autosome+".gz"), "w")

pops = numpy.array([pop.strip() for pop in inpop], dtype = "U100")
mafs = numpy.zeros(pops.shape, dtype = "U25")
numpy.savetxt(outfile, pops[None], delimiter = " ", fmt="%.21s")

c = 0
line = infile.readline()
oldrs = None
for n, line in enumerate(infile):
    line = line.strip().split()
    if autosome is not None and autosome != line[0]: continue
    rs = line[1]
    if oldrs is None:
        oldrs = line[1]
    if oldrs != rs:
        numpy.savetxt(outfile, mafs[None], delimiter = " ",  fmt="%.21s")
        mafs = numpy.zeros(pops.shape, dtype = "U25")
        oldrs = rs
        c += 1
    mafs[pops == line[2]] = ','.join([line[6], str(int(line[7])-int(line[6]))])

print "Processed {} SNPs".format(c)
        

