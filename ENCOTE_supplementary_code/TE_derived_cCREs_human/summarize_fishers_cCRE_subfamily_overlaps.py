'''
Gets the bedtools fisher's test results for all cCRE-TE subfamily overlaps in a given directory.
Multiple test corrects using Bonferonni method.
Usage: python3 summarize_fishers_cCRE_subfamily_overlaps.py <input directory of cCRE-subfamily tests> <output file>
'''

import sys, os

if len(sys.argv) != 3:
	sys.exit(__doc__)

dir_path = sys.argv[1].strip('/') + '/'

di = {}
list_overlaps = []
list_subfamilies = []
list_ccres = []
list_files = os.listdir(dir_path)
total_comparisons = len(list_files)
for filename in list_files:
	subfamily = '_'.join(filename.split('_')[:-1])
	ccre = filename.split('_')[-1].split('.')[0]
	if subfamily not in list_subfamilies:
		list_subfamilies.append(subfamily)
	if ccre not in list_ccres:
		list_ccres.append(ccre)
	with open(dir_path + filename, 'r') as f:
		for line in f:
			if 'db intervals' in line:
				total = line.split(':')[1].strip()
				continue
			elif 'overlaps' in line:
				num_overlaps = line.split(':')[1].strip()
				continue
			elif line.startswith('#'):
				continue
			elif line.startswith('left'): #p-value header line
				continue
			#At this point, it should be p-values from fishers
			fields = line.rstrip('\n').split()
			left = min([float(fields[0])*total_comparisons, 1])
			right = min([float(fields[1])*total_comparisons, 1])
			two_tailed = min([float(fields[2])*total_comparisons, 1])
			ratio = fields[3]
			overlap = subfamily + '_' + ccre
			di[overlap] = [subfamily, ccre, total, num_overlaps, str(left), str(right), str(two_tailed), ratio]
			list_overlaps.append(overlap)

with open(sys.argv[-1], 'w') as o:
	o.write('Subfamily\tcCRE-type\tSubfamily_number\tNum_overlaps\tLeft_padj\tRight_padj\t2tail_padj\tRatio\n')
	#Write corrected fishers tests to output file. Doing it this way to group same subfamily tests together
	for subfamily in list_subfamilies:
		for ccre in list_ccres:
			overlap = subfamily + '_' + ccre
			if overlap in list_overlaps:
				o.write('\t'.join(di[overlap]) + '\n')

