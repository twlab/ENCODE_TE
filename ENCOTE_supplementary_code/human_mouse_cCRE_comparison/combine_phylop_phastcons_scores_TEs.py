'''
Combines phyloP and phastCons scores for same regions from bedtools intersect.
Uses mean score (weighted by number of bp covered) to combine.
This version is for TE overlap.
Usage: python3 combine_phylop_phastcons_scores_TEs.py <intersect file> <output file>
'''

import sys

if len(sys.argv) != 3:
	sys.exit(__doc__)

elements = []
di = {}
scores = {}
with open(sys.argv[1], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		element = fields[0] + ':' + fields[1] + '-' + fields[2]
		element_info = fields[:6]
		if element not in di:
			elements.append(element)
			di[element] = element_info
		if fields[10] != '.':
			score = float(fields[10])
		else:
			print(line.rstrip('\n'))
			continue
		bp_overlap = int(fields[11])
		if element not in scores:
			scores[element] = [0, 0]
		scores[element][0] += score * bp_overlap
		scores[element][1] += bp_overlap

with open(sys.argv[-1], 'w') as o:
	for element in elements:
		new_line = '\t'.join(di[element])
		if element not in scores: #No phyloP/phastCons score for element
			new_line += '\tNA\t0\t0\n'
			o.write(new_line)
			continue
		element_len = int(di[element][2]) - int(di[element][1])
		mean_score = scores[element][0] / scores[element][1]
		coverage = scores[element][1] / element_len
		new_line += '\t' + str(mean_score) + '\t' + str(scores[element][1]) + '\t' + str(coverage) + '\n'
		o.write(new_line)

