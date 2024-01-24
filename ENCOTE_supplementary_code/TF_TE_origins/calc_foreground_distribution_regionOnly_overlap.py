'''
Calculates consensus coverage for cCRE overlapping (foreground) elements in general and the cCRE overlapping region only.
Outputs coverages for foreground/overlap regions to <output file basename> separated by TE subfamily.
Prints KS-test comparing full foreground coverage vs. cCRE overlapping regions only coverage to standard output.
This version only considers the region of the TE that overlaps a cCRE.
Usage: python3 calc_foreground_distribution_regionOnly_overlap.py <list of TE subfamilies> <directory of region overlap files> <directory of alignFastas> <output file basename>
'''

import sys, os
import numpy as np
from scipy.stats import ks_2samp

exclude_thresh = 2 #Threshold for minimum foreground/background enrichment in coverage required to exclude motifs in the region
exclude_base = 10 #Threshold for minimum foreground coverage required to exclude motifs in the region

if len(sys.argv) != 5:
	sys.exit(__doc__)

fdr = 0.05 #False-discovery rate for adjusted p-value threshold

def p_adjust_bh(p):
	p = np.asfarray(p)
	by_descend = p.argsort()[::-1]
	by_orig = by_descend.argsort()
	steps = float(len(p)) / np.arange(len(p), 0, -1)
	q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
	return q[by_orig]

list_subfamilies = []
with open(sys.argv[1], 'r') as f:
	for line in f:
		subfamily = line.rstrip('\n')
		list_subfamilies.append(subfamily)

overlap_dir = sys.argv[2].rstrip('/') + '/'
all_overlap_files = os.listdir(overlap_dir)
di_region = {}
for subfamily in list_subfamilies:
	region_file = subfamily + '.region_wInfo'
	if region_file in all_overlap_files:
		di_region[subfamily] = {}
		di_region[subfamily]['total'] = 0
		with open(overlap_dir + region_file, 'r') as f:
			for line in f:
				fields = line.rstrip('\n').split('\t')
				name = fields[0] + ':' + fields[1] + '-' + fields[2] + '(' + fields[5] + ')'
				strand = fields[5]
				if name not in di_region[subfamily]:
					di_region[subfamily][name] = []
				TE_start = int(fields[1])
				TE_stop = int(fields[2])
				ccre_start = int(fields[8])
				ccre_stop = int(fields[9])
				overlap_start = max([TE_start, ccre_start])
				overlap_stop = min([TE_stop, ccre_stop])
				if strand == '+': #Forward strand
					start = overlap_start - TE_start
					stop = overlap_stop - TE_start
				else: #Reverse strand
					start = abs(overlap_stop - TE_stop)
					stop = abs(overlap_start - TE_stop)
				di_region[subfamily][name].append([start, stop])
				di_region[subfamily]['total'] += 1

align_dir = sys.argv[3].rstrip('/') + '/'
list_align_files = os.listdir(align_dir)
full_coverage_di = {}
overlap_coverage_di = {}
for subfamily in list_subfamilies:
	all_full_coverage = []
	all_overlap_coverage = []
	filename = subfamily + '.alignFasta'
	if filename in list_align_files:
		if subfamily in di_region:
			with open(align_dir + filename, 'r') as f:
				for line in f:
					if line.startswith('>'):
						name = line.lstrip('>').rstrip('\n')
						seq = f.readline().rstrip('\n')
						ref = f.readline().rstrip('\n')
						if name not in di_region[subfamily]:
							continue
						sorted_overlaps = sorted(di_region[subfamily][name])
						overlap_coverage = []
						for i in range(len(sorted_overlaps)):
							overlap_coverage.append([])
						full_coverage = []
						seq_base = 0
						for i in range(len(seq)):
							if seq[i] == '-': #Deletion in element => no coverage at reference base
								full_coverage.append(0)
								for j in range(len(sorted_overlaps)):
									overlap_coverage[j].append(0)
							else:
								if ref[i] == '-': #Insertion in element => no reference base
									seq_base += 1
									continue
								else: #Base in element and reference
									full_coverage.append(1)
									for j in range(len(sorted_overlaps)):
										if seq_base >= sorted_overlaps[j][0] and seq_base <= sorted_overlaps[j][1]:
											overlap_coverage[j].append(1)
										else:
											overlap_coverage[j].append(0)
									seq_base += 1
						#Add coverage for each element and each overlapping region
						if np.any(all_full_coverage):
							all_full_coverage = np.add(all_full_coverage, full_coverage)
							for i in range(len(sorted_overlaps)):
								all_overlap_coverage = np.add(all_overlap_coverage, overlap_coverage[i])
						else:
							all_full_coverage = full_coverage
							all_overlap_coverage = overlap_coverage[0]
							for i in range(len(sorted_overlaps)-1):
								all_overlap_coverage = np.add(all_overlap_coverage, overlap_coverage[i+1])
				full_coverage_di[subfamily] = all_full_coverage
				overlap_coverage_di[subfamily] = all_overlap_coverage

list_tested_subfamilies = []
kstest_pval = []
kstest_stat = []
for subfamily in list_subfamilies:
	if subfamily in full_coverage_di:
		all_full_coverage = full_coverage_di[subfamily]
		all_overlap_coverage = overlap_coverage_di[subfamily]
		if all_full_coverage == []: #No foreground elements
			sys.stderr.write(subfamily + '\thas no foreground elements\n')
			continue
		with open(sys.argv[-1] + '_' + subfamily, 'w') as o:
			num_elements = len(di_region[subfamily]) - 1 #Subtract 1 for the 'total' key
			num_overlaps = di_region[subfamily]['total']
			o.write('Subfamily\tPosition\tOverlap\tForeground\tEnrichment\n')
			for i in range(len(all_overlap_coverage)):
				if all_full_coverage[i] == 0:
					enrichment = 'NA'
				else:
					enrichment = round((all_overlap_coverage[i]/num_overlaps)/(all_full_coverage[i]/num_elements), 4)
				o.write(subfamily + '\t' + str(i) + '\t' + str(all_overlap_coverage[i]) + '\t' + str(all_full_coverage[i]) + '\t' + str(enrichment) + '\n')
		#Do KS-test
		full_forKS = []
		overlap_forKS = []
		for i in range(len(all_overlap_coverage)):
			full_forKS.extend([i]*all_full_coverage[i])
			overlap_forKS.extend([i]*all_overlap_coverage[i])
		try:
			stat, pval = ks_2samp(overlap_forKS, full_forKS)
		except ValueError:
			print(subfamily)
			print(all_full_coverage)
			print(all_overlap_coverage)
			raise
		list_tested_subfamilies.append(subfamily)
		kstest_pval.append(pval)
		kstest_stat.append(stat)

kstest_qval = p_adjust_bh(kstest_pval)

print('Subfamily\tp-adj\tp-val\tstatistic')
for i in range(len(list_tested_subfamilies)):
	subfamily = list_tested_subfamilies[i]
	print(subfamily + '\t' + str(kstest_qval[i]) + '\t' + str(kstest_pval[i]) + '\t' + str(kstest_stat[i]))

