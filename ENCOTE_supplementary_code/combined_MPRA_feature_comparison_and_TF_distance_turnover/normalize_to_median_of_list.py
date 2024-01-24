'''
Normalize expression levels (columns 7 to end/9) by the median of a given list of controls.
Controls file must be element name and expression (tab-separated).
Usage: python3 normalize_to_median_of_list.py <input common file format> <controls file> <output file>
'''

import sys
import statistics as stats

if len(sys.argv) != 4:
	sys.exit(__doc__)

controls = []
with open(sys.argv[2], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		value = float(fields[1])
		controls.append(value)

norm = stats.median(controls)

with open(sys.argv[-1], 'w') as o:
	with open(sys.argv[1], 'r') as f:
		for line in f:
			fields = line.rstrip('\n').split('\t')
			new_line = '\t'.join(fields[:6])
			for exp in fields[6:]:
				new_line += '\t' + str(float(exp) - norm)
			o.write(new_line + '\n')

