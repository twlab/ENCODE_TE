'''
Summarizes consensus coverage info for significant TE subfamilies per TE family.
Takes a directory of coverage files as input.
Requires file with list of input TE subfamilies (subfamilies to include).
Usage: python3 summarize_consensus_coverage_inputSubfamilies_byFamily.py <input coverage directory> <list input TE subfamilies> <TE classes file> <output file basename>
'''

import sys, os
import statistics as stats
import math

if len(sys.argv) != 5:
	sys.exit(__doc__)

di_subfamilies = []
with open(sys.argv[2], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		di_subfamilies.append(subfamily)

list_families = []
TE_info = {'family': {}, 'class': {}}
with open(sys.argv[3], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		if subfamily not in di_subfamilies:
			continue
		family = fields[1]
		TE_info['family'][subfamily] = family
		if family not in list_families:
			list_families.append(family)
		if "DNA" in family:
			TEclass = "DNA"
		elif "LTR" in family:
			if "-int" in subfamily or "MamGypsy2" in subfamily:
				TEclass = "ERV-int"
			else:
				TEclass = "LTR"
		elif "SINE" in family:
			TEclass = "SINE"
		elif "LINE" in family:
			TEclass = "LINE"
		elif "SVA" in family:
			TEclass = "SVA"
		else:
			TEclass = "Other"
		TE_info['class'][subfamily] = TEclass

coverage_dir = sys.argv[1].rstrip('/') + '/'
try:
	list_files = os.scandir(path=coverage_dir)
except OSError:
	sys.stderr.write('Something wrong with opening <input coverage directory> ' + sys.argv[1])
	raise
	sys.exit()

di = {}
for filename in list_files:
	if not filename.is_file():
		continue
	name = filename.name
	if '.lengths' in name: #Ignore lengths info file
		continue
	subfamily = name.replace('cCRE_coverage_', '') #Get subfamily
	if subfamily in di:
		sys.stderr.write(subfamily + ' found twice in ' + sys.argv[1])
		sys.exit()
	di[subfamily] = {}
	with open(filename.path, 'r') as f:
		for line in f:
			fields = line.rstrip('\n').split('\t')
			annotation = fields[0]
			if annotation == 'enrichment':
				di[subfamily][annotation] = [float(x) for x in fields[1:]]
			else:
				di[subfamily][annotation] = [int(x) for x in fields[1:]]
	di[subfamily]['length'] = len(di[subfamily]['enrichment'])

di_results_family = {}
di_results_class = {}
#Get coverage info for all input TE subfamilies
for subfamily in sorted(di):
	if subfamily not in di_subfamilies:
		continue
	#Decide which base positions to consider for analysis based on number of elements in foreground and background that cover the base (don't want very low numbers)
	consensus_len = di[subfamily]['length']
	bases_to_use = []
	total = 0
	for i in range(consensus_len):
		if di[subfamily]['background'][i] >= 10: #Background elements have at least 10 coverage at the base
			bases_to_use.append(i)
			total += di[subfamily]['enrichment'][i]
	if len(bases_to_use) < 0.5*di[subfamily]['length']: #Less than 50% of bases have enough background coverage, skip subfamily and report
		sys.stderr.write(subfamily + '\thas too few useable bases (not enough background element coverage).\n')
		continue
	mean_enrichment = total/len(bases_to_use)
	#Condense enrichment across consensus to 100 bins (normalizes length to 100bp)
	bin_di = {}
	factor = consensus_len/100
	for base in bases_to_use:
		consensus_bin = math.floor(base/factor)
		if consensus_bin not in bin_di:
			bin_di[consensus_bin] = []
		bin_di[consensus_bin].append(di[subfamily]['enrichment'][base]/mean_enrichment) #Normalize enrichment to mean enrichment across consensus length of the TE subfamily
	#Add to results
	if subfamily not in TE_info['family']:
		TE_info['family'][subfamily] = 'Other'
		TE_info['class'][subfamily] = 'Other'
	family = TE_info['family'][subfamily]
	TEclass = TE_info['class'][subfamily]
	family_class = family + '_' + TEclass
	if family not in di_results_family:
		di_results_family[family] = []
		for consensus_bin in range(100):
			if consensus_bin in bin_di:
				di_results_family[family].append([stats.mean(bin_di[consensus_bin])])
			else:
				di_results_family[family].append([])
	else:
		for consensus_bin in range(100):
			if consensus_bin in bin_di:
				di_results_family[family][consensus_bin].append(stats.mean(bin_di[consensus_bin]))
	if family_class not in di_results_class:
		di_results_class[family_class] = []
		for consensus_bin in range(100):
			if consensus_bin in bin_di:
				di_results_class[family_class].append([stats.mean(bin_di[consensus_bin])])
			else:
				di_results_class[family_class].append([])
	else:
		for consensus_bin in range(100):
			if consensus_bin in bin_di:
				di_results_class[family_class][consensus_bin].append(stats.mean(bin_di[consensus_bin]))

with open(sys.argv[-1] + '_FamilyOnly', 'w') as o:
	o.write('Family\tConsensus_bin\tEnrichment\tUpper\tLower\n')
	for family in di_results_family:
		for consensus_bin in range(100):
			try:
				o.write(family + '\t' + str(consensus_bin) + '\t' + str(stats.mean(di_results_family[family][consensus_bin])) + '\t' +
					str(stats.mean(di_results_family[family][consensus_bin])+stats.stdev(di_results_family[family][consensus_bin])) + '\t' +
					str(stats.mean(di_results_family[family][consensus_bin])-stats.stdev(di_results_family[family][consensus_bin])) + '\n')
			except stats.StatisticsError: #No data in the bin across all TE subfamilies in the family
				o.write(family + '\t' + str(consensus_bin) + '\tNA\tNA\tNA\n')

with open(sys.argv[-1] + '_FamilyClass', 'w') as o:
	o.write('Family\tConsensus_bin\tEnrichment\tUpper\tLower\n')
	for family_class in di_results_class:
		for consensus_bin in range(100):
			try:
				o.write(family_class + '\t' + str(consensus_bin) + '\t' + str(stats.mean(di_results_class[family_class][consensus_bin])) + '\t' +
					str(stats.mean(di_results_class[family_class][consensus_bin])+stats.stdev(di_results_class[family_class][consensus_bin])) + '\t' +
					str(stats.mean(di_results_class[family_class][consensus_bin])-stats.stdev(di_results_class[family_class][consensus_bin])) + '\n')
			except stats.StatisticsError: #No data in the bin across all TE subfamilies in the family-class
				o.write(family_class + '\t' + str(consensus_bin) + '\tNA\tNA\tNA\n')

