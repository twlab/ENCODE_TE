'''
Summarizes phastCons scores for syntenic TF binding sites.
Usage: python3 calc_syntenicPeaks_phastCons_summary.py <overlap phastCons file> <output file>
'''

import sys
import numpy as np

if len(sys.argv) != 3:
	sys.exit(__doc__)

all_TFs = {}
ref_phastcons = {}
syn_phastcons = {}
with open(sys.argv[1], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		TF = fields[2]
		if TF not in all_TFs:
			all_TFs[TF] = 'Y'
		if TF not in ref_phastcons:
			ref_phastcons[TF] = []
			syn_phastcons[TF] = []
		if fields[-2] != 'NA':
			ref_phastcons[TF].append(float(fields[-2]))
		if fields[-1] != 'NA':
			syn_phastcons[TF].append(float(fields[-1]))

with open(sys.argv[-1], 'w') as o:
	o.write('TF\tGenome\tMean\tMedian\tMax\tPerc95\tNumber\n')
	for TF in sorted(all_TFs):
		if TF in ref_phastcons:
			o.write(TF + '\tReference\t' + str(np.mean(ref_phastcons[TF])) + '\t' + str(np.median(ref_phastcons[TF])) + '\t' + str(max(ref_phastcons[TF])) + '\t' + 
				str(np.quantile(ref_phastcons[TF], 0.95)) + '\t' + str(len(ref_phastcons[TF])) + '\n')
		else:
			o.write(TF + '\tReference\tNA\tNA\tNA\tNA\t0\n')
		if TF in syn_phastcons:
			o.write(TF + '\tSyntenic\t' + str(np.mean(syn_phastcons[TF])) + '\t' + str(np.median(syn_phastcons[TF])) + '\t' + str(max(syn_phastcons[TF])) + '\t' + 
				str(np.quantile(syn_phastcons[TF], 0.95)) + '\t' + str(len(syn_phastcons[TF])) + '\n')
		else:
			o.write(TF + '\tSyntenic\tNA\tNA\tNA\tNA\t0\n')

