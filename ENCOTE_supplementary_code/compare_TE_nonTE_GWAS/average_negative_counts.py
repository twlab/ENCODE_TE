'''
Calculate the average (mean) of negative control (shuffled) counts.
Usage: python3 average_negative_counts.py <counts files> <output file>
'''

import sys
import statistics as stats

if len(sys.argv) < 4:
	sys.exit(__doc__)

list_files = sys.argv[1:-1]
di = {}
for filename in list_files:
	with open(filename, 'r') as f:
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

with open(sys.argv[-1], 'w') as o:
	o.write(header)
	for parent_term in sorted(di):
		o.write(parent_term + '\t' + str(stats.mean(di[parent_term]['ccre'])) + '\t' + str(stats.mean(di[parent_term]['nonTE'])) + '\t' + str(stats.mean(di[parent_term]['TE'])) + '\t' +
			str(stats.mean(di[parent_term]['human_nonTE'])) + '\t' + str(stats.mean(di[parent_term]['human_TE'])) + '\n')

