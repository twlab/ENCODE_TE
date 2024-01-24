'''
Calculates summary statistics/numbers for variant frequency in shuffles.
	Mean, percent varFreq=0
Usage: python3 calc_summary_stats_shuffles.py <shuffle varFreq files> <output file basename>
'''

import sys
import statistics as stats

if len(sys.argv) < 4:
	sys.exit(__doc__)

list_files = sys.argv[1:-1]

shuffle_results_nonZero = {}
shuffle_results_includeZero = {}
for filename in list_files:
	di = {'TE': {'all': [], 'PLS': [], 'pELS': [], 'dELS': [], 'DNase-H3K4me3': [], 'CTCF-only': []},
		'non-TE': {'all': [], 'PLS': [], 'pELS': [], 'dELS': [], 'DNase-H3K4me3': [], 'CTCF-only': []}}
	with open(filename, 'r') as f:
		for line in f:
			fields = line.rstrip('\n').split('\t')
			varFreq = float(fields[-1])
			ccre = fields[5].replace(',CTCF-bound', '')
			if fields[6] == 'NA':
				di['non-TE']['all'].append(varFreq)
				di['non-TE'][ccre].append(varFreq)
			else:
				di['TE']['all'].append(varFreq)
				di['TE'][ccre].append(varFreq)
	shuffle_results_nonZero[filename] = {}
	shuffle_results_includeZero[filename] = {}
	for annotation in sorted(di):
		for ccre in sorted(di[annotation]):
			shuffle_results_includeZero[filename][annotation + ccre] = [annotation, ccre, str(stats.mean(di[annotation][ccre])), 
											str(di[annotation][ccre].count(0)/len(di[annotation][ccre]))]
			shuffle_results_nonZero[filename][annotation + ccre] = [annotation, ccre, str(stats.mean([x for x in di[annotation][ccre] if x != 0])), 
										str(di[annotation][ccre].count(0)/len(di[annotation][ccre]))]

with open(sys.argv[-1] + '_includeZero', 'w') as o:
	o.write('Shuffle\tAnnotation\tcCRE_type\tMean\tPerc_0\n')
	for filename in sorted(shuffle_results_includeZero):
		for annotation in sorted(shuffle_results_includeZero[filename]):
			o.write(filename + '\t' + '\t'.join(shuffle_results_includeZero[filename][annotation]) + '\n')

with open(sys.argv[-1] + '_nonZero', 'w') as o:
	o.write('Shuffle\tAnnotation\tcCRE_type\tMean\tPerc_0\n')
	for filename in sorted(shuffle_results_nonZero):
		for annotation in sorted(shuffle_results_nonZero[filename]):
			o.write(filename + '\t' + '\t'.join(shuffle_results_nonZero[filename][annotation]) + '\n')

