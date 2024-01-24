'''
Gets TE subfamilies with length bias in cCRE vs. non-cCRE elements for KS tests.
No correction for multiple hypothesis testing (conservative, rejects more subfamilies).
Prints TE subfamilies that do not have length bias to standard output.
Usage: python3 get_length_biased_subfamilies_kstest_noCorrect.py <KS-test results file> <output file>
'''

import sys
#import numpy as np

if len(sys.argv) != 3:
	sys.exit(__doc__)

fdr = 0.05 #False discovery rate, adjust if needed

#Not used
#def p_adjust_bh(p):
#	p = np.asfarray(p)
#	by_descend = p.argsort()[::-1]
#	by_orig = by_descend.argsort()
#	steps = float(len(p)) / np.arange(len(p), 0, -1)
#	q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
#	return q[by_orig]

order = []
di = {}
non = {}
with open(sys.argv[1], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		order.append(subfamily)
		if fields[1] == 'NA':
			continue
		pval = float(fields[1])
		di[subfamily] = pval

with open(sys.argv[-1], 'w') as o:
	for subfamily in order:
		if subfamily not in di:
			print(subfamily)
			continue
		pval = di[subfamily]
		if pval < fdr:
			o.write(subfamily + '\n')
		else:
			print(subfamily)

