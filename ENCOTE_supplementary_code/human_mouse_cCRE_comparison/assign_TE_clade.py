'''
Assign TE clade info to TE-cCREs.
Usage: python3 assign_TE_clade.py <input cCRE info file> <shared TE clades file> <lineage TE clades file> <output file>
'''

import sys

if len(sys.argv) != 5:
	sys.exit(__doc__)

di = {}
di['shared'] = {}
with open(sys.argv[2], 'r') as f:
	header = f.readline().rstrip('\n').split('\t')
	clade_index = header.index('Clade')
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		clade = fields[clade_index]
		di['shared'][subfamily] = clade
di['lineage'] = {}
with open(sys.argv[3], 'r') as f:
	header = f.readline().rstrip('\n').split('\t')
	clade_index = header.index('Clade')
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		clade = fields[clade_index]
		di['lineage'][subfamily] = clade

with open(sys.argv[1], 'r') as f:
	with open(sys.argv[-1], 'w') as o:
		header = f.readline().rstrip('\n')
		o.write(header + '\tSpecies_classification\tClade\n')
		for line in f:
			fields = line.rstrip('\n').split('\t')
			subfamily = fields[5]
			if subfamily == 'NA':
				continue
			if subfamily in di['shared']:
				o.write(line.rstrip('\n') + '\tShared_subfamily\t' + di['shared'][subfamily] + '\n')
			elif subfamily in di['lineage']:
				o.write(line.rstrip('\n') + '\tLineage_subfamily\t' + di['lineage'][subfamily] + '\n')
			else:
				o.write(line.rstrip('\n') + '\tNA\tNA\n')

