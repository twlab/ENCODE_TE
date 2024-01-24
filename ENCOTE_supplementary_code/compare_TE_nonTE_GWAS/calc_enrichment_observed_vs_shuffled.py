'''
Calculates enrichment of observed vs. shuffled for counts of GWAS variant overlap with cCREs.
Prints empirical p-value/results for observed counts based on shuffled counts (background/expectation) for each GWAS parent term.
Usage: python3 calc_enrichment_observed_vs_shuffled.py <observed counts file> <directory of shuffled counts files> <output file>
'''

import sys, os

if len(sys.argv) != 4:
	sys.exit(__doc__)

dir_shuffled = sys.argv[2].rstrip('/') + '/'
list_files = os.listdir(dir_shuffled)
di = {}
for filename in list_files:
	with open(dir_shuffled + filename, 'r') as f:
		header = f.readline()
		for line in f:
			fields = line.rstrip('\n').split('\t')
			parent_term = fields[0]
			if parent_term not in di:
				di[parent_term] = {'ccre': [], 'nonTE': [], 'TE': [], 'human_nonTE': [], 'human_TE': []}
			di[parent_term]['ccre'].append(int(fields[1]))
			di[parent_term]['nonTE'].append(int(fields[2]))
			di[parent_term]['TE'].append(int(fields[3]))
			di[parent_term]['human_nonTE'].append(int(fields[4]))
			di[parent_term]['human_TE'].append(int(fields[5]))

observed = {}
with open(sys.argv[1], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		parent_term = fields[0]
		observed[parent_term] = {'ccre': int(fields[1]), 'nonTE': int(fields[2]), 'TE': int(fields[3]), 'human_nonTE': int(fields[4]), 'human_TE': int(fields[5])}

with open(sys.argv[-1], 'w') as o:
	o.write('Parent_term\tcCRE_group\tEnrichment\tShuffled_count\n')
	for parent_term in sorted(observed):
		for group in sorted(observed[parent_term]):
			for i in range(len(di[parent_term][group])):
				o.write(parent_term + '\t' + group + '\t' + str(observed[parent_term][group]/di[parent_term][group][i]) + '\t' + str(di[parent_term][group][i]) + '\n')

print('Parent_term\tcCRE_group\tNumber\tp-value')
for parent_term in sorted(observed):
	for group in sorted(observed[parent_term]):
		number = 0
		observed_count = observed[parent_term][group]
		total = len(di[parent_term][group])
		for i in range(total):
			if observed_count <= di[parent_term][group][i]:
				number += 1
		print(parent_term + '\t' + group + '\t' + str(number) + '\t' + str(number/total))

