'''
Takes matched mouse cCRE to human cCRE/TE files to quantify syntenic elements.
Adds age of TE subfamily.
Usage: python3 quantify_mouse_human_ccre_TEs.py <matched cCRE file> <matched TE file> <kimura divergences mouse> <kimura divergences human> <output file>
'''

import sys

if len(sys.argv) != 6:
	sys.exit(__doc__)

mouse_ages = {}
with open(sys.argv[3], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		age = fields[-1]
		mouse_ages[subfamily] = age

human_ages = {}
with open(sys.argv[4], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		age = fields[-1]
		human_ages[subfamily] = age

list_mm10_ccres = []
di = {}
with open(sys.argv[1], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		mm10_ccre = fields[0]
		list_mm10_ccres.append(mm10_ccre)
		di[mm10_ccre] = fields[:9]
		if fields[4] != "NA":
			subfamily = fields[5]
			if subfamily in mouse_ages:
				di[mm10_ccre].append(mouse_ages[subfamily])
			else:
				di[mm10_ccre].append('NA')
		else:
			di[mm10_ccre].append('NA')
		di[mm10_ccre] += fields[9:]

with open(sys.argv[2], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		mm10_ccre = fields[0]
		di[mm10_ccre] += fields[10:]
		if fields[10] != "NA":
			subfamily = fields[11]
			if subfamily in human_ages:
				di[mm10_ccre].append(human_ages[subfamily])
			else:
				di[mm10_ccre].append('NA')
		else:
			di[mm10_ccre].append('NA')

mouse_cCREs = {}
mouse_TEs = {}
syntenic_regions = []
human_cCREs = {} #must be subset of num_syntenic regions
human_TEs = {} #must be subset of num_syntenic regions
human_ccre_types = {}
mouse_ccre_types = {}
for mm10_ccre in list_mm10_ccres:
	fields = di[mm10_ccre]
	#Mouse cCRE info
	mouse_ccre_type = []
	if 'dELS' in fields[3]:
		mouse_ccre_type.append('dELS')
	if 'pELS' in fields[3]:
		mouse_ccre_type.append('pELS')
	if 'PLS' in fields[3]:
		mouse_ccre_type.append('PLS')
	if 'H3K4me3-DNase' in fields[3]:
		mouse_ccre_type.append('H3K4me3-DNase')
	if 'CTCF-only' in fields[3]:
		mouse_ccre_type.append('CTCF-only')
	if 'CTCF' in fields[3]:
		mouse_ccre_type.append('CTCF-bound')
	if 'DNase-only' in fields[3]:
		mouse_ccre_type.append('DNase-only')
	for ccre_type in mouse_ccre_type:
		if ccre_type not in mouse_ccre_types:
			mouse_ccre_types[ccre_type] = []
		mouse_ccre_types[ccre_type].append(mm10_ccre)
	mouse_cCREs[mm10_ccre] = mouse_ccre_type
	#Mouse TE, if any
	if fields[4] != 'NA':
		subfamily = fields[5]
		TE_class = fields[6].split('/')[0]
		try:
			family = fields[6].split('/')[1]
		except IndexError:
			family = 'NA'
		mouse_TEs[mm10_ccre] = [subfamily, family, TE_class]
	#Syntenic human region, if any
	if fields[10] != 'NA':
		syntenic_regions.append(mm10_ccre)
	#Syntenic human cCRE, if any
	if fields[14] != 'NA':
		human_ccre_type = []
		if 'dELS' in fields[14]:
			human_ccre_type.append('dELS')
		if 'pELS' in fields[14]:
			human_ccre_type.append('pELS')
		if 'PLS' in fields[14]:
			human_ccre_type.append('PLS')
		if 'H3K4me3-DNase' in fields[14]:
			human_ccre_type.append('H3K4me3-DNase')
		if 'CTCF-only' in fields[14]:
			human_ccre_type.append('CTCF-only')
		if 'CTCF' in fields[14]:
			human_ccre_type.append('CTCF-bound')
		if 'DNase-only' in fields[14]:
			human_ccre_type.append('DNase-only')
		for ccre_type in human_ccre_type:
			if ccre_type not in human_ccre_types:
				human_ccre_types[ccre_type] = []
			human_ccre_types[ccre_type].append(mm10_ccre)
		human_cCREs[mm10_ccre] = human_ccre_type
	#Syntenic human TE, if any
	if fields[16] != 'NA':
		subfamily = fields[17]
		TE_class = fields[18].split('/')[0]
		try:
			family = fields[18].split('/')[1]
		except IndexError:
			family = 'NA'
		human_TEs[mm10_ccre] = [subfamily, family, TE_class]

#Output
with open(sys.argv[-1], 'w') as o:
	o.write('Number mouse cCREs: ' + str(len(mouse_cCREs)) + '\n')
	o.write('Number mouse TE-cCREs: ' + str(len(mouse_TEs)) + '\n')
	o.write('Number syntenic regions in human: ' + str(len(syntenic_regions)) + '\n')
	o.write('Number syntenic human cCREs: ' + str(len(human_cCREs)) + '\n')
	o.write('Number syntenic human TEs: ' + str(len(human_TEs)) + '\n')
	for mm10_ccre in list_mm10_ccres:
		o.write('\t'.join(di[mm10_ccre]) + '\n')

