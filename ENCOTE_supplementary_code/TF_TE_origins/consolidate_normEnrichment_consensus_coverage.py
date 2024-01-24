'''
Combines consensus coverage normalized enrichment across different randomizations.
Usage: python3 consolidate_normEnrichment_consensus_coverage.py <normalized enrichment files> <output file>
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
			if fields[2] == 'NA':
				continue
			TEclass = fields[0]
			if TEclass not in di:
				di[TEclass] = {}
			consensus_bin = fields[1]
			if consensus_bin not in di[TEclass]:
				di[TEclass][consensus_bin] = {'mean': [], 'upper': [], 'lower': []}
			di[TEclass][consensus_bin]['mean'].append(float(fields[2]))
			di[TEclass][consensus_bin]['upper'].append(float(fields[3]))
			di[TEclass][consensus_bin]['lower'].append(float(fields[4]))

with open(sys.argv[-1], 'w') as o:
	o.write(header)
	for TEclass in sorted(di):
		for consensus_bin in sorted(di[TEclass]):
			o.write(TEclass + '\t' + consensus_bin + '\t' + str(stats.mean(di[TEclass][consensus_bin]['mean'])) + '\t' +
				str(stats.mean(di[TEclass][consensus_bin]['upper'])) + '\t' + str(stats.mean(di[TEclass][consensus_bin]['lower'])) + '\n')

