'''
Gets the mean percentage of elements with motif for each TE subfamily.
Usage: python3 get_mean_percent_elements_wMotifs.py <input per Motif percentage file> <output file>
'''

import sys
import statistics as stats

if len(sys.argv) != 3:
	sys.exit(__doc__)

di = {}
info = {}
with open(sys.argv[1], 'r') as f:
	header = f.readline().rstrip('\n').split('\t')
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		archetype = fields[4]
		percent = float(fields[6])
		if subfamily not in di:
			di[subfamily] = {'percent': [], 'archetype': []}
			info[subfamily] = fields[:4]
		di[subfamily]['percent'].append(percent)
		di[subfamily]['archetype'].append(archetype)

with open(sys.argv[-1], 'w') as o:
	o.write('\t'.join(header[:4]) + '\tPercent_motifs\tMotif_clusters\n')
	for subfamily in sorted(di):
		o.write('\t'.join(info[subfamily]) + '\t' + str(stats.mean(di[subfamily]['percent'])) + '\t' + ','.join(di[subfamily]['archetype']) + '\n')

