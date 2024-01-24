'''
Finds the motifs that enriched in intersected elements relative to non-intersected elements.
Intersect and non-intersect files can be output from bedtools intersect (-v option for non-intersect)
Uses Fisher's exact test to determine enrichment. Multiple hypothesis correct p-values using BH method.
Contingency table:		+intersect	-intersect
	elements_wMotif		    ##		    ##
	elements_without	    ##		    ##
Usage: python3 find_motifs_assoc_wIntersect_fishers_v3.py <input stringent motifs file> <input motif annotation file> <input intersect file> <input non-intersect file> <output file>
'''


import sys
import numpy as np
import scipy.stats

if len(sys.argv) != 6:
	sys.exit(__doc__)

fdr = 0.05 #False discovery rate is 0.05, change as necessary

def p_adjust_bh(p):
	p = np.asfarray(p)
	by_descend = p.argsort()[::-1]
	by_orig = by_descend.argsort()
	steps = float(len(p)) / np.arange(len(p), 0, -1)
	q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
	return q[by_orig]

stringent_motifs = []
with open(sys.argv[1], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		motif = fields[2]
		stringent_motifs.append(motif)

di_motifs = {}
order_elements = []
with open(sys.argv[2], 'r') as f:
	header = f.readline().rstrip('\n')
	fields = header.split('\t')
	list_motifs = []
	for motif in fields[2:]:
		motif_name = motif
		list_motifs.append(motif_name)
	num_motifs = len(list_motifs)
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		name = fields[1]
		di_motifs[name] = {}
		for i in range(num_motifs):
			motif = list_motifs[i]
			di_motifs[name][motif] = fields[2+i]
		order_elements.append(name)

di_fishers_cat = {}
with open(sys.argv[3], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		name = fields[0] + ':' + fields[1] + '-' + fields[2] + '(' + fields[5] + ')'
		di_fishers_cat[name] = {}
		element_motifs = di_motifs[name]
		if element_motifs == {}: #No motifs scanned for this subfamily, skip
			continue
		for motif in element_motifs:
			if int(element_motifs[motif]) >= 1: #Has at least 1 motif
				di_fishers_cat[name][motif] = '+motif+annotation'
			else: #Does not have motif
				di_fishers_cat[name][motif] = '-motif+annotation'

with open(sys.argv[4], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		name = fields[0] + ':' + fields[1] + '-' + fields[2] + '(' + fields[5] + ')'
		if name in di_fishers_cat: #TE in both intersect and non-intersect files, should not happen
			print(name + ' in intersect file ' + sys.argv[2] + ' and non-intersect file ' + sys.argv[3])
			print('Elements in the two files must be mutually exclusive')
			sys.exit(__doc__)
		di_fishers_cat[name] = {}
		element_motifs = di_motifs[name]
		if element_motifs == {}: #No motifs scanned for this subfamily, skip
			continue
		for motif in element_motifs:
			if int(element_motifs[motif]) >= 1: #Has at least 1 motif
				di_fishers_cat[name][motif] = '+motif-annotation'
			else: #Does not have motif
				di_fishers_cat[name][motif] = '-motif-annotation'	

with open(sys.argv[-1], 'w') as o:
	o.write('Subfamily\tMotif\tOdds_ratio\tadj-pval\tp-value\t+Motif+Annotation\t+Motif-Annotation\t-Motif+Annotation\t-Motif-Annotation\n')
	fishers_tests = []
	list_pvals = []
	category_counts = {}
	for motif in stringent_motifs:
		category_counts[motif] = [0, 0, 0, 0]
	#Get counts for contingency tables for each motif
	for name in order_elements:
		if name not in di_fishers_cat: #Skip if element not in either intersect nor non-intersect files (can happen if background element not selected)
			continue
		for motif in stringent_motifs:
			try:
				category = di_fishers_cat[name][motif]
			except KeyError:
				print(subfamily + ' does not have annotation for ' + motif + '\tSkipping...')
				stringent_motifs.remove(motif)
				continue
			if category == '+motif+annotation':
				category_counts[motif][0] += 1
			elif category == '+motif-annotation':
				category_counts[motif][1] += 1
			elif category == '-motif+annotation':
				category_counts[motif][2] += 1
			elif category == '-motif-annotation':
				category_counts[motif][3] += 1
	#Perform fisher's exact tests for each motif
	for motif in stringent_motifs:
		odds, pval = scipy.stats.fisher_exact([category_counts[motif][:2], category_counts[motif][2:]], alternative='greater')
		list_pvals.append(pval)
		fishers_tests.append([motif, str(odds), str(category_counts[motif][0]), str(category_counts[motif][1]), str(category_counts[motif][2]), str(category_counts[motif][3])])
	#Multiple hypothesis correct using Benjamini-Hochberg method
	list_padj = p_adjust_bh(list_pvals)
	for i in range(len(fishers_tests)):
		fishers_test = fishers_tests[i]
		padj = list_padj[i]
		new_line = subfamily + '\t' + '\t'.join(fishers_test[:2]) + '\t' + str(padj) + '\t' + str(list_pvals[i]) + '\t' + '\t'.join(fishers_test[2:]) + '\n'
		o.write(new_line)

