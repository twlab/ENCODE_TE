'''
Consolidates cCRE overlap region only coverage (relative to full TE sequences that overlap cCRE).
Usage: python3 consolidate_ccreOnly_overlap_coverage.py <directory of coverage file> <list of TE subfamilies file> <TE classes file> <output file>
'''

import sys, os
import math
import statistics as stats

if len(sys.argv) != 5:
	sys.exit(__doc__)

di_subfamilies = {}
all_families = {}
with open(sys.argv[3], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		family = fields[1]
		di_subfamilies[subfamily] = family
		if family not in all_families:
			all_families[family] = 'Y'
list_subfamilies = []
with open(sys.argv[2], 'r') as f:
	for line in f:
		subfamily = line.rstrip('\n')
		if subfamily not in di_subfamilies:
			di_subfamilies[subfamily] = 'Other'
			all_families['Other'] = 'Y'
		list_subfamilies.append(subfamily)

coverage_dir = sys.argv[1].rstrip('/') + '/'
all_files = os.listdir(coverage_dir)
di_coverage = {}
for subfamily in list_subfamilies:
	coverage_file = 'coverage_' + subfamily
	if coverage_file in all_files:
		di_coverage[subfamily] = {'foreground': [], 'enrichment': []}
		with open(coverage_dir + coverage_file, 'r') as f:
			header = f.readline()
			for line in f:
				fields = line.rstrip('\n').split('\t')
				position = fields[1]
				if fields[-1] == 'NA':
					enrichment = 'NA'
				else:
					enrichment = float(fields[-1])
				foreground_num = int(fields[3])
				di_coverage[subfamily]['enrichment'].append(enrichment)
				di_coverage[subfamily]['foreground'].append(foreground_num)

di_results = {}
#Get coverage info for all input TE subfamilies
for subfamily in list_subfamilies:
	if subfamily not in di_subfamilies or subfamily not in di_coverage:
		continue
	#Decide which base positions to consider for analysis based on number of elements in foreground and background that cover the base (don't want very low numbers)
	consensus_len = len(di_coverage[subfamily]['enrichment'])
	bases_to_use = []
	total = 0
	for i in range(consensus_len):
		if di_coverage[subfamily]['foreground'][i] >= 10: #Foreground elements have at least 10 coverage at the base
			bases_to_use.append(i)
			total += di_coverage[subfamily]['enrichment'][i]
	try:
		mean_enrichment = total/len(bases_to_use)
	except ZeroDivisionError:
		continue
	#Condense enrichment across consensus to 100 bins (normalizes length to 100bp)
	bin_di = {}
	factor = consensus_len/100
	for base in bases_to_use:
		consensus_bin = math.floor(base/factor)
		if consensus_bin not in bin_di:
			bin_di[consensus_bin] = []
		bin_di[consensus_bin].append(di_coverage[subfamily]['enrichment'][base]/mean_enrichment) #Normalize enrichment to mean enrichment across consensus length of the TE subfamily
	#Add to results
	family = di_subfamilies[subfamily]
	if family not in di_results:
		di_results[family] = []
		for consensus_bin in range(100):
			if consensus_bin in bin_di:
				di_results[family].append([stats.mean(bin_di[consensus_bin])])
			else:
				di_results[family].append([])
	else:
		for consensus_bin in range(100):
			if consensus_bin in bin_di:
				di_results[family][consensus_bin].append(stats.mean(bin_di[consensus_bin]))

with open(sys.argv[-1], 'w') as o:
	o.write('Family\tPosition\tEnrichment\n')
	for family in sorted(all_families):
		if family in di_results:
			for consensus_bin in range(100):
				try:
					o.write(family + '\t' + str(consensus_bin) + '\t' + str(stats.mean(di_results[family][consensus_bin])) + '\n')
				except stats.StatisticsError: #No data in the bin across all non-uniform TE subfamilies in the class
					o.write(family + '\t' + str(consensus_bin) + '\tNA\n')

