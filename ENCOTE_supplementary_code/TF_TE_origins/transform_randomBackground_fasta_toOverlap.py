'''
Transforms sequences picked for random background to overlap format (for coverage analysis)
Usage: python3 transform_randomBackground_fasta_toOverlap.py <input file> <output file>
'''

import sys

if len(sys.argv) != 3:
	sys.exit(__doc__)

with open(sys.argv[1], 'r') as f:
	with open(sys.argv[-1], 'w') as o:
		for line in f:
			if line.startswith('>'):
				name = line.rstrip('\n').lstrip('>')
				chrom  = name.split(':')[0]
				start = name.split(':')[1].split('-')[0]
				stop = name.split('-')[1].split('(')[0]
				strand = name.split('(')[1].rstrip(')')
				o.write(chrom + '\t' + start + '\t' + stop + '\tNA\tNA\t' + strand + '\n')

