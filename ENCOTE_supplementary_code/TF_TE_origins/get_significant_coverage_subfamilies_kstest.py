'''
Uses KS-test to test if foreground sequences have the same length distribution as background sequences.
	Specifically, the test is whether the foreground and background sequence lengths were drawn from the same distribution.
Prints KS-test result to standard output. Uses <group name> to specify what the KS-test result is for.
This version keeps all background elements.
Usage: python3 get_significant_coverage_subfamilies_kstest.py <input foreground fasta sequences> <input background fasta sequences> <group name>
'''

import sys
import math
from scipy.stats import ks_2samp

if len(sys.argv) != 4:
	sys.exit(__doc__)

foreground_di = {}
foreground = []
#min_len = 1000000
#max_len = 0
#foreground_total_len = 0
with open(sys.argv[1], 'r') as f:
	for line in f:
		if line.startswith('>'):
			name = line.rstrip('\n').lstrip('>')
			seq = f.readline().rstrip('\n')
			length = len(seq)
			if name in foreground_di:
				sys.stderr.write(name + ' found twice in ' + sys.argv[1] + '\tExiting...\n')
				sys.exit()
			foreground_di[name] = length
			foreground.append(length)
#			foreground_total_len += length
#			if length < min_len:
#				min_len = length
#			if length > max_len:
#				max_len = length
if len(foreground_di) == 0: #No foreground elements
	sys.stderr.write('No foreground elements for ' + sys.argv[1] + '\n')
	sys.exit()

background_di = {}
background = []
#background_total_len = 0
with open(sys.argv[2], 'r') as f:
	for line in f:
		if line.startswith('>'):
			name = line.rstrip('\n').lstrip('>')
			seq = f.readline().rstrip('\n')
			length = len(seq)
#			if length < min_len:
#				min_len = length
#			if length > max_len:
#				max_len = length
			if name in background_di:
				sys.stderr.write(name + ' found twice in ' + sys.argv[2] + '\tExiting...\n')
				sys.exit()
			background_di[name] = length
			background.append(length)
#			background_total_len += length

#Perform KS-test to determine whether length distribution is different between foreground and background
if len(background) == 0:
	sys.stderr.write('No background elements for ' + sys.argv[2] + '\n')
	sys.exit()
stat, pval = ks_2samp(foreground, background, alternative='less')
print(sys.argv[3] + '\t' + str(pval) + '\t' + str(stat))

