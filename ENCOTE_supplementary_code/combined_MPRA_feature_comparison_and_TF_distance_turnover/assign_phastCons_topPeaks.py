'''
Assigns phastCons scores to top peaks peaks.
Usage: python3 assign_phastCons_topPeaks.py <top peaks file> <overlaps with phastCons ref species> <overlaps with phastCons syntenic species> <output file>
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
		header = f.readline()
		o.write(header.rstrip() + '\tRef_phastCons\tSyntenic_phastCons\n')
		for line in f:
			fields = line.rstrip('\n').split('\t')
			ref_coord = fields[4].split('|')[0]
			syn_coords = [x.split('|')[0] for x in fields[5].split(',')]
			o.write(line.rstrip())
			new_line = ''
			if ref_coord in ref_di:
				new_line += '\t' + ref_di[ref_coord] + '\t'
			else:
				new_line += '\tNA\t'
			for syn_coord in syn_coords:
				if syn_coord in syn_di:
					new_line += syn_di[syn_coord] + ','
				else:
					new_line += 'NA,'
			o.write(new_line.rstrip(',') + '\n')

