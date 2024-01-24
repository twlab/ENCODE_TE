'''
Divides a Needle alignment file of multiple sequences into separate files of one alignment each.
Usage: python3 divide_needleAlign_individualFiles.py <input Needle outfile> <output basename>
'''

import sys

if len(sys.argv) != 3:
	sys.exit(__doc__)

count = 1
basename = sys.argv[2]
with open(sys.argv[1], 'r') as f:
	num_alignInfo_brackets = 0
	export_lines = ''
	for line in f:
		if '#===' in line:
			num_alignInfo_brackets += 1
			if num_alignInfo_brackets > 2:
				with open(basename + str(count), 'w') as o:
					o.write(export_lines)
				count += 1
				num_alignInfo_brackets -= 2
				export_lines = ''
		export_lines += line

#Export last alignment
with open(basename + str(count), 'w') as o:
	o.write(export_lines)

