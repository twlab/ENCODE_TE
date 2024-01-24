'''
Removes duplicate cCRE-TE overlaps due to repeatmasker having two TEs in the same genomic locus.
Usage: python3 remove_duplicate_overlaps.py <input file> <output file>
'''

import sys

if len(sys.argv) != 3:
	sys.exit(__doc__)

di = {}
ccres = []
with open(sys.argv[1], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		ccre_id = fields[3]
		if ccre_id not in ccres:
			ccres.append(ccre_id)
			di[ccre_id] = fields
		else:
			overlap = int(fields[-1])
			prev_overlap = int(di[ccre_id][-1])
			if overlap > prev_overlap:
				di[ccre_id] = fields

with open(sys.argv[-1], 'w') as o:
	for ccre_id in ccres:
		o.write('\t'.join(di[ccre_id]) + '\n')

