'''
Match TE subfamily info from tsv files.
First column of each tsv file must be the TE subfamily.
Usage: python3 match_subfamily_info.py <file 1> <file 2> <output file>
'''

import sys

if len(sys.argv) != 4:
	sys.exit(__doc__)

di = {}
subfamilies = []
with open(sys.argv[1], 'r') as f:
	line = f.readline()
	fields = line.rstrip('\n').split('\t')
	if fields[0].lower() == 'subfamily': #Header, first column is 'subfamily'
		header = fields
	else: #First line is not header
		subfamily = fields[0]
		di[subfamily] = fields[1:]
		subfamilies.append(subfamily)
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		di[subfamily] = fields[1:]
		subfamilies.append(subfamily)

shared = []
with open(sys.argv[2], 'r') as f:
	line = f.readline()
	fields = line.rstrip('\n').split('\t')
	if fields[0].lower() == 'subfamily': #Header, first column is 'subfamily'
		header += fields[1:]
	else: #First line is not header
		subfamily = fields[0]
		di[subfamily] += fields[1:]
		if subfamily in subfamilies:
			shared.append(subfamily)
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		if subfamily in di:
			di[subfamily] += fields[1:]
		else:
			continue
		if subfamily in subfamilies:
			shared.append(subfamily)

with open(sys.argv[-1], 'w') as o:
	o.write('\t'.join(header) + '\n')
	for subfamily in shared:
		o.write(subfamily + '\t' + '\t'.join(di[subfamily]) + '\n')

