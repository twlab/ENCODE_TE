'''
Calculates class-level log2 enrichments at element (count) and base pair levels.
	Features: Removes separation of CTCF-bound/unbound cCREs (e.g. "dELS,CTCF-bound" and "dELS" are combined)
		  Focuses on 4 main human classes (DNA, LTR, LINE, SINE) + SVA. Combines subcategories into those main groups + an "Other" group
Print TE class overlap numbers to standard output.
Usage: python3 get_class_enrichments_cCRE_overlap_combined.py <TE rmsk bed file> <input overlap file> <output file>
'''

import sys
import math

if len(sys.argv) != 4:
	sys.exit(__doc__)

genome_size = 3209286105 #~3.2 billion bp genome size, number from adding up chromosome sizes from /bar/genomes/hg38/hg38.chrom.sortV.sizes

di_TEs = {}
list_TE_names = []
with open(sys.argv[1], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		start = int(fields[1])
		stop = int(fields[2])
		TE_length = stop - start
		if "DNA" in fields[4]:
			TE_name = "DNA"
		elif "LTR" in fields[4]:
			TE_name = "LTR"
		elif "LINE" in fields[4]:
			TE_name = "LINE"
		elif "SINE" in fields[4]:
			TE_name = "SINE"
		elif "SVA" in fields[4]:
			TE_name = "SVA"
		else:
			TE_name = "Other"
		if TE_name in di_TEs:
			di_TEs[TE_name][0] += 1
			di_TEs[TE_name][1] += TE_length
		else:
			di_TEs[TE_name] = [1, TE_length]
			list_TE_names.append(TE_name)

ccre_di = {}
ccre_types = []
overlap_di = {}
with open(sys.argv[2], 'r') as f:
	header_fields = f.readline().rstrip('\n').split('\t')
	ccre_index = 'NA'
	name_index = 'NA'
	overlap_index = 'NA'
	for i in range(len(header_fields)):
		if header_fields[i] == 'cCRE-type':
			ccre_index = i
		elif header_fields[i] == 'repClass/Family':
			name_index = i
		elif header_fields[i] == 'bpOverlap':
			overlap_index = i
	if ccre_index == 'NA' or name_index == 'NA' or overlap_index == 'NA':
		print(sys.argv[2] + ' does not have header with cCRE-type, repName, bpOverlap fields')
		sys.exit(__doc__)
	prev_TE = ''
	for line in f:
		fields = line.rstrip('\n').split('\t')
		start = int(fields[1])
		stop = int(fields[2])
		ccre_length = stop - start
		ccre_type = fields[ccre_index].split(',')[0]
		if ccre_type in ccre_di:
			ccre_di[ccre_type][0] += 1
			ccre_di[ccre_type][1] += ccre_length
		else:
			ccre_di[ccre_type] = [1, ccre_length]
			ccre_types.append(ccre_type)
		if fields[name_index] == '.': #No TE overlap, go to next line
			continue
		elif "DNA" in fields[name_index]:
			TE_name = "DNA"
		elif "LTR" in fields[name_index]:
			TE_name = "LTR"
		elif "LINE" in fields[name_index]:
			TE_name = "LINE"
		elif "SINE" in fields[name_index]:
			TE_name = "SINE"
		elif "SVA" in fields[name_index]:
			TE_name = "SVA"
		else:
			TE_name = "Other"
		bpOverlap = int(fields[overlap_index])
		index = ccre_types.index(ccre_type)
		if TE_name not in overlap_di:
			overlap_di[TE_name] = []
		while len(overlap_di[TE_name]) < len(ccre_types):
			overlap_di[TE_name].append([0, 0])
		TE = fields[6] + ':' + fields[7] + '-' + fields[8]
		if TE != prev_TE: #Not the same TE as the previous TE, add 1 to cCRE overlap count
			overlap_di[TE_name][index][0] += 1
		overlap_di[TE_name][index][1] += bpOverlap
		prev_TE = TE

with open(sys.argv[-1], 'w') as o:
	#Header
	o.write('Class')
	for ccre_type in ccre_types:
		o.write('\t' + ccre_type + '(count)\t' + ccre_type + '(bases)')
	o.write('\n')
	#Class enrichments
	for TE_name in list_TE_names:
		#No overlap of class elements with any cCRE
		if TE_name not in overlap_di:
			o.write(TE_name)
			for ccre_type in ccre_types:
				o.write('\tNA\tNA')
			o.write('\n')
			continue
		#Some overlap of class elements with at least one cCRE
		o.write(TE_name)
		TE_overlap_info = overlap_di[TE_name]
		while len(TE_overlap_info) < len(ccre_types):
			TE_overlap_info.append([0, 0])
		for i in range(len(TE_overlap_info)):
			TE_ccre_info = TE_overlap_info[i]
			ccre_type = ccre_types[i]
			if TE_ccre_info[0] == 0: #No overlap of the TE class with the cCRE
				o.write('\tNA\tNA')
			else:
				count_enrich = math.log2((float(TE_ccre_info[0])/di_TEs[TE_name][1])/(float(ccre_di[ccre_type][0])/genome_size))
				base_enrich = math.log2((float(TE_ccre_info[1])/di_TEs[TE_name][1])/(float(ccre_di[ccre_type][1])/genome_size))
				o.write('\t' + str(count_enrich) + '\t' + str(base_enrich))
		o.write('\n')
	#All Classes excluding Other (DNA, LTR, LINE, SINE, SVA) enrichments
	all_TE_ccres = []
	while len(all_TE_ccres) < len(ccre_types):
		all_TE_ccres.append([0, 0])
	all_TE_genome = [0, 0]
	for TE_name in list_TE_names:
		if TE_name == "Other":
			continue
		all_TE_genome[0] += di_TEs[TE_name][0]
		all_TE_genome[1] += di_TEs[TE_name][1]
		if TE_name not in overlap_di:
			continue
		TE_overlap_info = overlap_di[TE_name]
		while len(TE_overlap_info) < len(ccre_types):
			TE_overlap_info.append([0, 0])
		for i in range(len(TE_overlap_info)):
			all_TE_ccres[i][0] += TE_overlap_info[i][0]
			all_TE_ccres[i][1] += TE_overlap_info[i][1]
	o.write('TE')
	for i in range(len(all_TE_ccres)):
		TE_ccre_info = all_TE_ccres[i]
		ccre_type = ccre_types[i]
		if TE_ccre_info[0] == 0: #No overlap of any TE with the cCRE
			o.write('\tNA\tNA')
		else:
			count_enrich = math.log2((float(TE_ccre_info[0])/all_TE_genome[1])/(float(ccre_di[ccre_type][0])/genome_size))
			base_enrich = math.log2((float(TE_ccre_info[1])/all_TE_genome[1])/(float(ccre_di[ccre_type][1])/genome_size))
			o.write('\t' + str(count_enrich) + '\t' + str(base_enrich))
	o.write('\n')

#Overlap info to standard output
header = 'Class\tClass_count\tClass_bases'
for ccre_type in ccre_types:
	header += '\t' + ccre_type + '(count)\t' + ccre_type + '(bases)'
print(header)
#Class counts
for TE_name in list_TE_names:
	line = TE_name + '\t' + str(di_TEs[TE_name][0]) + '\t' + str(di_TEs[TE_name][1])
	#No overlap of class elements with any cCRE
	if TE_name not in overlap_di:
		for ccre_type in ccre_types:
			line += '\t0\t0'
		print(line)
		continue
	#Some overlap of class elements with at least one cCRE
	TE_overlap_info = overlap_di[TE_name]
	while len(TE_overlap_info) < len(ccre_types):
		TE_overlap_info.append([0, 0])
	for TE_ccre_info in TE_overlap_info:
		line += '\t' + str(TE_ccre_info[0]) + '\t' + str(TE_ccre_info[1])
	print(line)

