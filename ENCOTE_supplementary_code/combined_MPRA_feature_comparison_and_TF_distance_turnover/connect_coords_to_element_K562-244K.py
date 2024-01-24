'''
Connects coordinates to expression levels for the K562 lentiMPRA dataset from the Ahituv lab.
Prints elements missing coordinates to standard output.
Usage: python3 connect_coords_to_element_K562_lentiMPRA.py <coordinates file> <expression file> <output file>
'''

import sys

if len(sys.argv) != 4:
	sys.exit(__doc__)

di_coords = {}
with open(sys.argv[1], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		name = fields[3]
		di_coords[name] = fields

with open(sys.argv[-1], 'w') as o:
	with open(sys.argv[2], 'r') as f:
		for line in f:
			fields = line.rstrip('\n').split('\t')
			name = fields[0]
			if name in di_coords: #should always be the case
				new_line = '\t'.join(di_coords[name]) + '\t' + fields[1]
				o.write(new_line + '\n')
			else:
				print(name)


