'''
Assigns TE subfamily age to file.
Usage: python3 assign_age.py <input file> <age file> <output file>
'''

import sys

if len(sys.argv) != 4:
	sys.exit(__doc__)

ages = {}
with open(sys.argv[2], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		age = fields[2]
		ages[subfamily] = age

with open(sys.argv[1], 'r') as f:
	with open(sys.argv[-1], 'w') as o:
		header = f.readline().rstrip('\n')
		o.write(header + '\tAge\n')
		for line in f:
			fields = line.rstrip('\n').split('\t')
			subfamily = fields[0]
			if subfamily in ages:
				o.write(line.rstrip('\n') + '\t' + ages[subfamily] + '\n')
			else:
				o.write(line.rstrip('\n') + '\tNA\n')

