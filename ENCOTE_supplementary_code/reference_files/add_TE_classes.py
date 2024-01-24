'''
Adds TE class info to subfamily clade/age info.
Usage: python3 add_TE_classes.py <subfamily clade/age file> <TE class file> <output file>
'''

import sys

if len(sys.argv) != 4:
	sys.exit(__doc__)

di = {}
with open(sys.argv[2], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		TE_class = fields[1]
		di[subfamily] = TE_class

with open(sys.argv[1], 'r') as f:
	with open(sys.argv[-1], 'w') as o:
		header = f.readline().rstrip('\n').split('\t')
		o.write(header[0] + '\tClass\t' + '\t'.join(header[1:]) + '\n')
		for line in f:
			fields = line.rstrip('\n').split('\t')
			subfamily = fields[0]
			if subfamily not in di:
				continue
			TE_class = di[subfamily]
			o.write(subfamily + '\t' + TE_class + '\t' + '\t'.join(fields[1:]) + '\n')

