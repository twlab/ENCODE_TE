'''
Assigns phastCons scores to syntenic peaks.
Usage: python3 assign_phastCons_syntenicPeaks.py <syntenic peaks file> <overlaps with phastCons ref species> <overlaps with phastCons syntenic species> <output file>
'''

import sys

if len(sys.argv) != 5:
	sys.exit(__doc__)

ref_di = {}
with open(sys.argv[2], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		coord = fields[0] + ':' + fields[1] + '-' + fields[2]
		phastcons = fields[3]
		ref_di[coord] = phastcons
syn_di = {}
with open(sys.argv[3], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		coord = fields[0] + ':' + fields[1] + '-' + fields[2]
		phastcons = fields[3]
		syn_di[coord] = phastcons

with open(sys.argv[1], 'r') as f:
	with open(sys.argv[-1], 'w') as o:
		o.write('Ref_peak\tSyntenic_region\tTF\tRef_phastCons\tSyntenic_phastCons\n')
		for line in f:
			fields = line.rstrip('\n').split('\t')
			syntenic_peak = fields[1].split('|')[-1]
			if syntenic_peak != 'yesPeak':
				continue
			ref_coord = fields[0].split('|')[0]
			syn_coord = fields[1].split('|')[0]
			o.write(line.rstrip() + '\t' + ref_di[ref_coord] + '\t' + syn_di[syn_coord] + '\n')

