'''
Get the number of cCRE overlaps for each element in a subfamily.
Usage: python3 get_ccre_overlaps_per_element_v2.py <directory of overlap files> <output file>
'''

import sys
import os

if len(sys.argv) != 3:
	sys.exit(__doc__)

dir_path = sys.argv[1].rstrip('/') + '/'
list_files = os.listdir(dir_path)
di_ccre = {}
di_info = {}
list_elements = []
list_ccres = []
for filename in list_files:
	with open(dir_path + filename, 'r') as f:
		celltype = '_'.join('.'.join(filename.split('.')[:-1]).split('_')[:-1])
		ccre = '.'.join(filename.split('.')[:-1]).split('_')[-1]
		if ccre not in list_ccres:
			list_ccres.append(ccre)
		ccre_index = list_ccres.index(ccre)
		subfamily = filename.split('.')[-1]
		for line in f:
			fields = line.rstrip('\n').split('\t')
			element = fields[0] + ':' + fields[1] + '-' + fields[2]
			if element not in list_elements:
				di_info[element] = fields[:6]
				di_ccre[element] = []
				while len(di_ccre[element]) < len(list_ccres):
					di_ccre[element].append(0)
				list_elements.append(element)
			if fields[7] != '.': #Overlap with a cCRE
				while len(di_ccre[element]) < len(list_ccres):
					di_ccre[element].append(0)
				di_ccre[element][ccre_index] += 1

with open(sys.argv[-1], 'w') as o:
	o.write('chrom\tstart\tstop\tSubfamily\tClass\tstrand')
	for ccre in list_ccres:
		o.write('\t' + ccre)
	o.write('\n')
	for element in list_elements:
		while len(di_ccre[element]) < len(list_ccres):
			di_ccre[element].append(0)
		o.write('\t'.join(di_info[element]) + '\t' + '\t'.join([str(x) for x in di_ccre[element]]) + '\n')

