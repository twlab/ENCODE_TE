'''
Assign phylop and phastcons scores to cCREs.
Usage: python3 assign_phylop_phastcons_scores_ccres.py <cCRE annotation file> <phylop scores file> <phastcons scores file> <output file>
'''

import sys

if len(sys.argv) != 5:
	sys.exit(__doc__)

di = {}
with open(sys.argv[2], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		coord = fields[0] + ':' + fields[1] + '-' + fields[2]
		id1 = fields[3]
		di[coord] = {}
		di[coord]['id'] = id1
		di[coord]['phylop'] = fields[6] + '\t' + fields[7] + '\t' + fields[8]
with open(sys.argv[3], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		coord = fields[0] + ':' + fields[1] + '-' + fields[2]
		id1 = fields[3]
		if coord not in di:
			di[coord] = {}
			di[coord]['id'] = id1
		else:
			if di[coord]['id'] != id1: #Check cCRE ID match
				sys.stderr.write('<phylop score file> and <phastcons scores file> unexpectedly have different IDs for ' + coord + '\nExiting...')
				sys.exit()
		di[coord]['phastcons'] = fields[6] + '\t' + fields[7] + '\t' + fields[8]

with open(sys.argv[1], 'r') as f:
	with open(sys.argv[-1], 'w') as o:
		header = f.readline().rstrip('\n')
		o.write(header + '\tphyloP_score\tphyloP_bases\tphyloP_coverage\tphastCons_score\tphastCons_bases\tphastCons_coverage\n')
		for line in f:
			o.write(line.rstrip('\n'))
			fields = line.rstrip('\n').split('\t')
			coord = fields[0]
			id1 = fields[1]
			if coord in di:
				if di[coord]['id'] != id1: #Check cCRE ID match
					sys.stderr.write('<cCRE annotation file> unexpectedly has different ID for ' + coord + '\nExiting...')
					sys.exit()
				if 'phylop' in di[coord]:
					o.write('\t' + di[coord]['phylop'])
				else:
					o.write('\tNA\tNA\tNA')
				if 'phastcons' in di[coord]:
					o.write('\t' + di[coord]['phastcons'] + '\n')
				else:
					o.write('\tNA\tNA\tNA\n')

