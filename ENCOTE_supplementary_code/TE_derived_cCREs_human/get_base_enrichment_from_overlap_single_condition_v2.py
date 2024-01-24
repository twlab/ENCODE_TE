'''
Gets the enrichment for each overlap file in a given directory.
Assumes each overlap file has the full TE subfamily (bedtools intersect -wao option with TE subfamily bed as -a file) and a single cCRE is overlapped.
Usage: python3 get_base_enrichment_from_overlap_single_condition_v2.py <input overlap directory> <all cCRE bed files directory> <output file>
'''

import sys, os
import math
import gzip

genome_size = 3209286105 #~3.2 billion bp genome size, number from adding up chromosome sizes from /bar/genomes/hg38/hg38.chrom.sortV.sizes

if len(sys.argv) != 4:
	sys.exit(__doc__)

ccre_path = sys.argv[2].rstrip('/') + '/'
list_ccre_files = os.listdir(ccre_path)
ccre_di = {}
for filename in list_ccre_files: #Assuming cCRE bed files in directory are gzipped (being lazy)
	celltype = '.'.join(filename.split('.')[:-2])
	with gzip.open(ccre_path + filename, mode='rt') as f:
		for line in f:
			ctcf = False
			fields = line.rstrip('\n').split('\t')
			if 'dELS' in fields[-2]:
				ccre = 'dELS'
			elif 'pELS' in fields[-2]:
				ccre = 'pELS'
			elif 'PLS' in fields[-2]:
				ccre = 'PLS'
			elif 'DNase-H3K4me3' in fields[-2]:
				ccre = 'DNase-H3K4me3'
			elif 'CTCF-only' in fields[-2]:
				ccre = 'CTCF-only'
			elif 'DNase-only' in fields[-2]:
				ccre = 'DNase-only'
			else: #Low-DNase, skip
				continue
			if 'CTCF' in fields[-2]:
				ctcf = True
			start = int(fields[1])
			stop = int(fields[2])
			ccre_len = stop - start
			key = celltype + '/' + ccre
			if key in ccre_di:
				ccre_di[key] += ccre_len
			else:
				ccre_di[key] = ccre_len
			if ctcf:
				key = celltype + '/CTCF-bound'
				if key in ccre_di:
					ccre_di[key] += ccre_len
				else:
					ccre_di[key] = ccre_len

dir_path = sys.argv[1].rstrip('/') + '/'
list_files = os.listdir(dir_path)
di = {}
list_keys = []
for filename in list_files:
	with open(dir_path + filename, 'r') as f:
		celltype = '_'.join('.'.join(filename.split('.')[:-1]).split('_')[:-1])
		ccre = '.'.join(filename.split('.')[:-1]).split('_')[-1]
		subfamily = filename.split('.')[-1]
		total_subfamily_len = 0
		overlap_len = 0
		num_overlap = 0
		num_elements = 0
		elements = [] #Keeps track of elements that overlap with cCRE, used so that elements are not double counted for overlap number
		for line in f:
			fields = line.rstrip('\n').split('\t')
			start = int(fields[1])
			stop = int(fields[2])
			te_len = stop - start
			total_subfamily_len += te_len
			num_elements += 1
			element = fields[0] + ':' + fields[1] + '-' + fields[2]
			if fields[7] != '.': #Overlap with a cCRE
				overlap_len += int(fields[-1])
				if element not in elements:
					num_overlap += 1
					elements.append(element)
		try:
			di[celltype + '/' + ccre] = [subfamily, celltype, ccre, str(num_elements), str(num_overlap), str(math.log2((overlap_len/total_subfamily_len)/(ccre_di[celltype + '/' + ccre]/genome_size)))]
		except ValueError:
			di[celltype + '/' + ccre] = [subfamily, celltype, ccre, str(num_elements), str(num_overlap), 'NA']
		list_keys.append(celltype + '/' + ccre)

with open(sys.argv[-1], 'w') as o:
	o.write('Subfamily\tCell/Tissue\tcCRE\tSubfamily_size\tnum_overlaps\tlog2_enrichment\n')
	for key in list_keys:
		o.write('\t'.join(di[key]) + '\n')

