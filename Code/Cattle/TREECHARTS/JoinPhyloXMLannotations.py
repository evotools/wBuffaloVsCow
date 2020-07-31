import sys
from colormap import rgb2hex
from Bio import Phylo


tree1 = Phylo.read(sys.argv[1],'phyloxml')
tree2 = Phylo.read(sys.argv[2],'phyloxml')
treef = tree1

n = 0
colors = []
m = 0
lengths = []

for clade in tree2.find_clades():
    m += 1
    lengths.append(clade.branch_length)

for n, clade in enumerate(treef.find_clades()):
    clade.branch_length = lengths[n]
Phylo.write(tree1, 'final.xml', 'phyloxml')
