import sys

outgroup = None
if len(sys.argv) == 3:
	if raw_input('Do you want to specify an outgroup (y/N): ').lower() == 'y':
		outgroup = raw_input('Specify outgroup (CASE SENSITIVE!): ')
		outgroup=outgroup.replace('_','').replace('-','')
elif len(sys.argv) == 4:
	outgroup = sys.argv[3]
	outgroup=outgroup.replace('_','').replace('-','')
else:
	sys.exit('Wrong command line:python PRG.py inputmatrix outgroup')

id_in = open(sys.argv[1] + '.id')
mat_in = open(sys.argv[1])



# Input collections
ids = []
fids = []
indexes = []
n_fids = {}
conv_fids = {}
conv_ids = {}
fid_idx = {}
nfid = 0
nind = 0

# Read id matrix
print 'Reading ID file...'
for n, i in enumerate(id_in):
	fid, iid = i.strip().split()
	ids.append(iid)
	fid = fid.replace('_', '').replace('-', '')
	fids.append(fid)	
	if outgroup != None and fid == outgroup:
		indexes.append(n + 1)
	if fid not in conv_fids:
		conv_fids[fid] = '%03i' % nfid
		fid_idx[fid] = 0
		nfid += 1
	if iid not in conv_ids:
		conv_ids[iid] = '%03i' % fid_idx[fid]
		fid_idx[fid] += 1
	nind += 1

print 'Done.\n Read %i populations.\n' % nfid

print 'Creating new ids for individuals...'
new_fids = [fid.replace('_', '')[0] + fid.replace('_', '')[2:4] if len(fid.split('_')[0]) > 3 else fid.replace('_','')[0:3] for fid in fids]
out_ids = [new_fids[e] + conv_fids[fids[e]] + conv_ids[ids[e]] for e in xrange(0, len(fids))]

# Create progressive new ids

out_conv = open('ID_convertion_file.txt', 'w')
out_mat = open('infile_{}', 'w')

print 'Writing convertion dataset in: ID_convertion_file.txt'
for n, i in enumerate(out_ids):
	out_conv.write('%s\t%s\t%s\n' % (fids[n], ids[n], i))
out_conv.close()

if outgroup != None:
	print 'Saving indexes of outgroup.'
	outgr = open('Outgroups_indexes.txt', 'w')
	[outgr.write('%i\n' % (i)) for i in indexes]
	outgr.close()


print 'Saving new matrix in: infile'
out_mat.write('%i\n' % nind)
for n, i in enumerate(mat_in):
	out_mat.write('%s %s' % (out_ids[n], i))
print 'Done.\n\nAll convertion completed.'

out_mat.close()



