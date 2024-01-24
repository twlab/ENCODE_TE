'''
Matches elements to the original input bed file after liftOver and bedtools intersect (-wao option to retain info about all liftOver-able elements).
Requires input bed file, liftOver output, liftOver unmapped file, and bedtools intersect file.
Usage: python3 match_elements_after_liftOver_intersect_hg38mm10.py <input bed file> <liftOver output file> <liftOver unmapped file> <bedtools intersect output file> <final output file>
'''

import sys

if len(sys.argv) != 6:
	sys.exit(__doc__)

#Get all input names before liftOver
in_order = []
in_di = {}
with open(sys.argv[1], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		name = fields[0] + ':' + fields[1] + '-' + fields[2]
		in_order.append(name)
		in_di[name] = line.rstrip('\n')

#Get input names that are unmapped during liftOver
unmapped = {}
with open(sys.argv[3], 'r') as f:
	for line in f:
		if line.startswith('#'):
			fields = f.readline().rstrip('\n').split()
			name = fields[0] + ':' + fields[1] + '-' + fields[2]
			unmapped[name] = 'Y'

#Get the association of name post-liftOver (lift_name) with the input name pre-liftOver (stored in liftOver_di)
mapped = {}
liftOver_di = {} #Used to direct the liftOver mapped name to the input name (pre-liftOver)
with open(sys.argv[2], 'r') as f:
	count_in = 0 #Keeps track of the input names in order
	for line in f:
		fields = line.rstrip('\n').split('\t')
		lift_name = fields[0] + ':' + fields[1] + '-' + fields[2]
		for i in range(count_in, len(in_order)):
			in_name = in_order[i]
			if in_name in unmapped: #if the input name was unmapped, go to the next input name
				continue
			else: #if the input name was mapped, make the association and update count_in
				if lift_name in mapped: #if input name maps to more than one lift_name
					liftOver_di[lift_name].append(in_name)
				else:
					liftOver_di[lift_name] = [in_name]
					mapped[lift_name] = 'Y'
				count_in = i+1
				break

added = {}
#Get the intersect info and associate back to the input names
with open(sys.argv[4], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		name = fields[0] + ':' + fields[1] + '-' + fields[2]
		if name in added:
			continue
		if name not in liftOver_di:
			continue
		if len(liftOver_di[name]) > 1: #More than one input name maps to the liftOver name
			for i in range(len(liftOver_di[name])):
				in_name = liftOver_di[name][i]
				in_di[in_name] = in_di[in_name] + '\t' + line.rstrip('\n')
		else: #One to one match
			in_name = liftOver_di[name][0]
			in_di[in_name] = in_di[in_name] + '\t' + line.rstrip('\n')
		added[name] = 'Y'

#Fill in empty spaces for unmapped input names using info from the intersect file
num_col = len(fields)
empty_line = '.\t-1\t-1'
for i in range(num_col-3):
	empty_line += '\t.'
for in_name in unmapped:
	in_di[in_name] = in_di[in_name] + '\t' + empty_line

#Write to final output file
with open(sys.argv[-1], 'w') as o:
	for name in in_order:
		o.write(in_di[name] + '\n')

