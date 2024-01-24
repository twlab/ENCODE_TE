'''
Get clades info from Dfam embl file.
Usage: python3 get_clades_dfam.py <Dfam embl file> <output file>
'''

import sys

if len(sys.argv) != 3:
	sys.exit(__doc__)

new_line = ''
with open(sys.argv[1], 'r') as f:
	with open(sys.argv[-1], 'w') as o:
		for line in f:
			if line.startswith('NM'):
				if new_line != '':
					o.write(new_line + '\n')
					new_line = ''
				new_line += line.lstrip('NM').strip()
			elif line.startswith('OS'):
				new_line += '\t' + line.lstrip('OS').strip()
		o.write(new_line + '\n')

