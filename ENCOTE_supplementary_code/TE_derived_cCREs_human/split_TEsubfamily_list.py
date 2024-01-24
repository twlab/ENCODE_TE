'''
Splits the list of TE subfamilies with X number of subfamilies per file, starting from a given TE subfamily.
Usage: python3 split_TEsubfamily_list.py <input file> <number of subfamilies per file> <starting TE subfamily> <output file basename>
'''

import sys
import math

if len(sys.argv) != 5:
	sys.exit(__doc__)

try:
	num_subfamilies = int(sys.argv[2])
except TypeError:
	print('<number of subfamilies per file> must be an integer')
	sys.exit(__doc__)

start = sys.argv[3]

subfamilies = []
with open(sys.argv[1], 'r') as f:
	for line in f:
		subfamilies.append(line.rstrip('\n'))

try:
	start_index = subfamilies.index(start)
except ValueError:
	print(start + ' not in list of TE subfamilies. Will begin from start of TE subfamilies file.')
	start_index = 0

num_files = math.ceil((len(subfamilies) - start_index) / num_subfamilies)

for i in range(num_files):
	filename = sys.argv[-1] + '_group' + str(i+1)
	with open(filename, 'w') as o:
		for j in range(num_subfamilies):
			index = start_index + i*num_subfamilies + j
			try:
				o.write(subfamilies[index] + '\n')
			except IndexError:
				sys.exit()

