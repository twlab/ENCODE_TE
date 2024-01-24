'''
This script is designed for comparing TE to non-TE cCREs and each to shuffled genomic regions, based on differences with flanks.
Calculates empirical p-values for observed variant frequencies based on shuffle summary statistics.
	1. Difference between TE and non-TE means vs. shuffle
	2. Difference between TE and non-TE percentage overlapping variants vs. shuffle
	3. Difference between TE/non-TE and flank mean vs. shuffle
	4. Difference between TE/non-TE and flank percentage overlapping variants vs. shuffle
This version excludes variant frequencies = 0 for mean calculation.
Usage: python3 calc_empirical_pvals_flanks_nonZero.py <observed varFreq file> <observed left flank varFreq file> <observed right flank varFreq file> 
							<shuffle summary stats file> <shuffle left flank summary stats file> <shuffle right flank summary stats file> <output file basename>
'''

import sys
import statistics as stats

if len(sys.argv) != 8:
	sys.exit(__doc__)

def get_observed_varFreq(filename):
	observed = {'TE': {'all': [], 'PLS': [], 'pELS': [], 'dELS': [], 'DNase-H3K4me3': [], 'CTCF-only': []},
			'non-TE': {'all': [], 'PLS': [], 'pELS': [], 'dELS': [], 'DNase-H3K4me3': [], 'CTCF-only': []}}
	numbers = {}
	with open(filename, 'r') as f:
		for line in f:
			fields = line.rstrip('\n').split('\t')
			ccre = fields[5].replace(',CTCF-bound', '')
			varFreq = float(fields[-1])
			if fields[6] == 'NA':
				observed['non-TE']['all'].append(varFreq)
				observed['non-TE'][ccre].append(varFreq)
			else:
				observed['TE']['all'].append(varFreq)
				observed['TE'][ccre].append(varFreq)
	for annotation in observed:
		numbers[annotation] = {}
		for ccre in observed[annotation]:
			numbers[annotation][ccre] = {}
			numbers[annotation][ccre]['mean'] = stats.mean([x for x in observed[annotation][ccre] if x != 0])
			numbers[annotation][ccre]['perc'] = 1 - observed[annotation][ccre].count(0)/len(observed[annotation][ccre])
	return numbers

def get_shuffled_summaries(filename):
	shuffled = {'TE': {'all': {}, 'PLS': {}, 'pELS': {}, 'dELS': {}, 'DNase-H3K4me3': {}, 'CTCF-only': {}},
			'non-TE': {'all': {}, 'PLS': {}, 'pELS': {}, 'dELS': {}, 'DNase-H3K4me3': {}, 'CTCF-only': {}}}
	for annotation in shuffled:
		for ccre in shuffled[annotation]:
			shuffled[annotation][ccre]['mean'] = {}
			shuffled[annotation][ccre]['perc'] = {}
	with open(filename, 'r') as f:
		header = f.readline()
		for line in f:
			fields = line.rstrip('\n').split('\t')
			shuffle_num = int(fields[0].replace('_noCDS', '').split('p')[-1]) #Shuffle number follows "overlap" after removing "_noCDS" from the end
			annotation = fields[1]
			ccre = fields[2]
			mean = float(fields[3])
			perc = 1 - float(fields[4])
			shuffled[annotation][ccre]['mean'][shuffle_num] = mean
			shuffled[annotation][ccre]['perc'][shuffle_num] = perc
	return shuffled

#Observed cCREs
observed_numbers = get_observed_varFreq(sys.argv[1])

#Observed cCRE flanks
observed_left = get_observed_varFreq(sys.argv[2])
observed_right = get_observed_varFreq(sys.argv[3])

#Shuffled cCRE coordinates
shuffled_stats = get_shuffled_summaries(sys.argv[4])
shuffled_left_stats = get_shuffled_summaries(sys.argv[5])
shuffled_right_stats = get_shuffled_summaries(sys.argv[6])

#Compare observed to shuffle
results = {}
list_ccres = sorted(observed_numbers['TE'])
#Difference between TE and non-TE means vs. shuffle
results['mean_diff'] = {}
with open(sys.argv[-1] + '_shuffled_mean_diff', 'w') as o:
	o.write('Shuffle_num\tcCRE_type\tMean_diff\n')
	for ccre in list_ccres:
		observed_diff = abs(observed_numbers['TE'][ccre]['mean'] - observed_numbers['non-TE'][ccre]['mean'])
		list_shuffled_diff = []
		for i in sorted(shuffled_stats['TE'][ccre]['mean']):
			shuffled_diff = abs(shuffled_stats['TE'][ccre]['mean'][i] - shuffled_stats['non-TE'][ccre]['mean'][i])
			list_shuffled_diff.append(shuffled_diff)
			o.write(str(i) + '\t' + ccre + '\t' + str(shuffled_diff) + '\n')
		counter = 0
		for i in range(len(list_shuffled_diff)):
			if observed_diff <  list_shuffled_diff[i]:
				counter += 1
		pval = counter/len(list_shuffled_diff)
		results['mean_diff'][ccre] = [pval, observed_diff, stats.mean(list_shuffled_diff)]
#Difference between TE and non-TE percentage of zeroes vs. shuffle
results['perc_diff'] = {}
with open(sys.argv[-1] + '_shuffled_perc_diff', 'w') as o:
	o.write('Shuffle_num\tcCRE_type\tPerc_diff\n')
	for ccre in list_ccres:
		observed_diff = abs(observed_numbers['TE'][ccre]['perc'] - observed_numbers['non-TE'][ccre]['perc'])
		list_shuffled_diff = []
		for i in sorted(shuffled_stats['TE'][ccre]['perc']):
			shuffled_diff = abs(shuffled_stats['TE'][ccre]['perc'][i] - shuffled_stats['non-TE'][ccre]['perc'][i])
			list_shuffled_diff.append(shuffled_diff)
			o.write(str(i) + '\t' + ccre + '\t' + str(shuffled_diff) + '\n')
		counter = 0
		for i in range(len(list_shuffled_diff)):
			if observed_diff <  list_shuffled_diff[i]:
				counter += 1
		pval = counter/len(list_shuffled_diff)
		results['perc_diff'][ccre] = [pval, observed_diff, stats.mean(list_shuffled_diff)]

#Difference between TE/non-TE and flank mean vs. shuffle
results['TE_mean_left'] = {}
results['TE_mean_right'] = {}
with open(sys.argv[-1] + '_shuffled_TE_mean', 'w') as o:
	o.write('Shuffle_num\tShuffle_dir\tcCRE_type\tDifference\n')
	for ccre in list_ccres:
		observed_mean_diff_left = abs(observed_numbers['TE'][ccre]['mean'] - observed_left['TE'][ccre]['mean'])
		observed_mean_diff_right = abs(observed_numbers['TE'][ccre]['mean'] - observed_right['TE'][ccre]['mean'])
		list_shuffled_diff_left = []
		list_shuffled_diff_right = []
		for i in sorted(shuffled_stats['TE'][ccre]['mean']):
			shuffled_diff_left = abs(shuffled_stats['TE'][ccre]['mean'][i] - shuffled_left_stats['TE'][ccre]['mean'][i])
			shuffled_diff_right = abs(shuffled_stats['TE'][ccre]['mean'][i] - shuffled_right_stats['TE'][ccre]['mean'][i])
			list_shuffled_diff_left.append(shuffled_diff_left)
			list_shuffled_diff_right.append(shuffled_diff_right)
			o.write(str(i) + '\tleft\t' + ccre + '\t' + str(shuffled_diff_left) + '\n')
			o.write(str(i) + '\tright\t' + ccre + '\t' + str(shuffled_diff_right) + '\n')
		counter_left = 0
		for i in range(len(list_shuffled_diff_left)):
			if observed_mean_diff_left < list_shuffled_diff_left[i]:
				counter_left += 1
		pval_left = counter_left/len(list_shuffled_diff_left)
		results['TE_mean_left'][ccre] = [pval_left, observed_mean_diff_left, stats.mean(list_shuffled_diff_left)]
		counter_right = 0
		for i in range(len(list_shuffled_diff_right)):
			if observed_mean_diff_right < list_shuffled_diff_right[i]:
				counter_right += 1
		pval_right = counter_right/len(list_shuffled_diff_right)
		results['TE_mean_right'][ccre] = [pval_right, observed_mean_diff_right, stats.mean(list_shuffled_diff_right)]
results['nonTE_mean_left'] = {}
results['nonTE_mean_right'] = {}
with open(sys.argv[-1] + '_shuffled_nonTE_mean', 'w') as o:
	o.write('Shuffle_num\tShuffle_dir\tcCRE_type\tDifference\n')
	for ccre in list_ccres:
		observed_mean_diff_left = abs(observed_numbers['non-TE'][ccre]['mean'] - observed_left['non-TE'][ccre]['mean'])
		observed_mean_diff_right = abs(observed_numbers['non-TE'][ccre]['mean'] - observed_right['non-TE'][ccre]['mean'])
		list_shuffled_diff_left = []
		list_shuffled_diff_right = []
		for i in sorted(shuffled_stats['non-TE'][ccre]['mean']):
			shuffled_diff_left = abs(shuffled_stats['non-TE'][ccre]['mean'][i] - shuffled_left_stats['non-TE'][ccre]['mean'][i])
			shuffled_diff_right = abs(shuffled_stats['non-TE'][ccre]['mean'][i] - shuffled_right_stats['non-TE'][ccre]['mean'][i])
			list_shuffled_diff_left.append(shuffled_diff_left)
			list_shuffled_diff_right.append(shuffled_diff_right)
			o.write(str(i) + '\tleft\t' + ccre + '\t' + str(shuffled_diff_left) + '\n')
			o.write(str(i) + '\tright\t' + ccre + '\t' + str(shuffled_diff_right) + '\n')
		counter_left = 0
		for i in range(len(list_shuffled_diff_left)):
			if observed_mean_diff_left < list_shuffled_diff_left[i]:
				counter_left += 1
		pval_left = counter_left/len(list_shuffled_diff_left)
		results['nonTE_mean_left'][ccre] = [pval_left, observed_mean_diff_left, stats.mean(list_shuffled_diff_left)]
		counter_right = 0
		for i in range(len(list_shuffled_diff_right)):
			if observed_mean_diff_right < list_shuffled_diff_right[i]:
				counter_right += 1
		pval_right = counter_right/len(list_shuffled_diff_right)
		results['nonTE_mean_right'][ccre] = [pval_right, observed_mean_diff_right, stats.mean(list_shuffled_diff_right)]

#Difference between TE/non-TE and flank percentage overlapping variants vs. shuffle
results['TE_perc_left'] = {}
results['TE_perc_right'] = {}
with open(sys.argv[-1] + '_shuffled_TE_perc', 'w') as o:
	o.write('Shuffle_num\tShuffle_dir\tcCRE_type\tDifference\n')
	for ccre in list_ccres:
		observed_perc_diff_left = abs(observed_numbers['TE'][ccre]['perc'] - observed_left['TE'][ccre]['perc'])
		observed_perc_diff_right = abs(observed_numbers['TE'][ccre]['perc'] - observed_right['TE'][ccre]['perc'])
		list_shuffled_diff_left = []
		list_shuffled_diff_right = []
		for i in sorted(shuffled_stats['TE'][ccre]['perc']):
			shuffled_diff_left = abs(shuffled_stats['TE'][ccre]['perc'][i] - shuffled_left_stats['TE'][ccre]['perc'][i])
			shuffled_diff_right = abs(shuffled_stats['TE'][ccre]['perc'][i] - shuffled_right_stats['TE'][ccre]['perc'][i])
			list_shuffled_diff_left.append(shuffled_diff_left)
			list_shuffled_diff_right.append(shuffled_diff_right)
			o.write(str(i) + '\tleft\t' + ccre + '\t' + str(shuffled_diff_left) + '\n')
			o.write(str(i) + '\tright\t' + ccre + '\t' + str(shuffled_diff_right) + '\n')
		counter_left = 0
		for i in range(len(list_shuffled_diff_left)):
			if observed_perc_diff_left < list_shuffled_diff_left[i]:
				counter_left += 1
		pval_left = counter_left/len(list_shuffled_diff_left)
		results['TE_perc_left'][ccre] = [pval_left, observed_perc_diff_left, stats.mean(list_shuffled_diff_left)]
		counter_right = 0
		for i in range(len(list_shuffled_diff_right)):
			if observed_perc_diff_right < list_shuffled_diff_right[i]:
				counter_right += 1
		pval_right = counter_right/len(list_shuffled_diff_right)
		results['TE_perc_right'][ccre] = [pval_right, observed_perc_diff_right, stats.mean(list_shuffled_diff_right)]
results['nonTE_perc_left'] = {}
results['nonTE_perc_right'] = {}
with open(sys.argv[-1] + '_shuffled_nonTE_perc', 'w') as o:
	o.write('Shuffle_num\tShuffle_dir\tcCRE_type\tDifference\n')
	for ccre in list_ccres:
		observed_perc_diff_left = abs(observed_numbers['non-TE'][ccre]['perc'] - observed_left['non-TE'][ccre]['perc'])
		observed_perc_diff_right = abs(observed_numbers['non-TE'][ccre]['perc'] - observed_right['non-TE'][ccre]['perc'])
		list_shuffled_diff_left = []
		list_shuffled_diff_right = []
		for i in sorted(shuffled_stats['non-TE'][ccre]['perc']):
			shuffled_diff_left = abs(shuffled_stats['non-TE'][ccre]['perc'][i] - shuffled_left_stats['non-TE'][ccre]['perc'][i])
			shuffled_diff_right = abs(shuffled_stats['non-TE'][ccre]['perc'][i] - shuffled_right_stats['non-TE'][ccre]['perc'][i])
			list_shuffled_diff_left.append(shuffled_diff_left)
			list_shuffled_diff_right.append(shuffled_diff_right)
			o.write(str(i) + '\tleft\t' + ccre + '\t' + str(shuffled_diff_left) + '\n')
			o.write(str(i) + '\tright\t' + ccre + '\t' + str(shuffled_diff_right) + '\n')
		counter_left = 0
		for i in range(len(list_shuffled_diff_left)):
			if observed_perc_diff_left < list_shuffled_diff_left[i]:
				counter_left += 1
		pval_left = counter_left/len(list_shuffled_diff_left)
		results['nonTE_perc_left'][ccre] = [pval_left, observed_perc_diff_left, stats.mean(list_shuffled_diff_left)]
		counter_right = 0
		for i in range(len(list_shuffled_diff_right)):
			if observed_perc_diff_right < list_shuffled_diff_right[i]:
				counter_right += 1
		pval_right = counter_right/len(list_shuffled_diff_right)
		results['nonTE_perc_right'][ccre] = [pval_right, observed_perc_diff_right, stats.mean(list_shuffled_diff_right)]

with open(sys.argv[-1] + '_results', 'w') as o:
	o.write('Comparison\tcCRE_type\tpval\tObserved\tShuffled\n')
	for comparison in sorted(results):
		for ccre in sorted(results[comparison]):
			o.write(comparison + '\t' + ccre + '\t' + '\t'.join([str(x) for x in results[comparison][ccre]]) + '\n')

