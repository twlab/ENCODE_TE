'''
Removes cCREs that have any overlap with coding sequence.
Usage: python3 remove_CDS_overlap_cCREs.py <cCRE annotation file> <cCRE CDS overlap> <output file>
'''

import sys

if len(sys.argv) != 4:
	sys.exit(__doc__)

di = {}
with open(sys.argv[2], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		ccre_id = fields[3]
		di[ccre_id] = 'Y'

with open(sys.argv[1], 'r') as f:
	with open(sys.argv[-1], 'w') as o:
		for line in f:
			fields = line.rstrip('\n').split('\t')
			ccre_id = fields[3]
			if ccre_id in di:
				continue
			o.write(line)

