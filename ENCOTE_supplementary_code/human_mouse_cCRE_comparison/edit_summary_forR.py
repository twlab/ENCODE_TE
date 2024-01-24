'''
Edits summary file for R plotting (because I hate R data formatting).
Usage: python3 edit_summary_forR.py <input file> <output file>
'''

import sys

if len(sys.argv) != 3:
	sys.exit(__doc__)

#Order for comparisons
order = ['same', 'different', 'human_only', 'mouse_only']

di = {}
with open(sys.argv[1], 'r') as f:
	for line in f:
		if line.startswith('#'):
			continue
		if line.rstrip('\n') == '':
			continue
		fields = line.rstrip('\n').split('\t')
		if fields[0] == 'Annotation':
			header = line
			continue
		annotation = fields[0]
		comparison = fields[1]
		ccre = fields[2]
		percent = float(fields[3])
		number = fields[4]
		if annotation not in di:
			di[annotation] = {}
		if ccre not in di[annotation]:
			di[annotation][ccre] = {}
		if comparison not in di[annotation][ccre]:
			di[annotation][ccre][comparison] = {}
		di[annotation][ccre][comparison]['percent'] = percent
		di[annotation][ccre][comparison]['number'] = number

for annotation in di:
	for ccre in di[annotation]:
		cumulative = 0
		for comparison in order:
			if comparison not in di[annotation][ccre]:
				continue
			percent = di[annotation][ccre][comparison]['percent']
			position = cumulative + percent/2
			di[annotation][ccre][comparison]['position'] = position
			cumulative += percent

with open(sys.argv[-1], 'w') as o:
	o.write(header.rstrip('\n') + '\tPosition\n')
	for annotation in di:
		for ccre in di[annotation]:
			for comparison in di[annotation][ccre]:
				o.write(annotation + '\t' + comparison + '\t' + ccre + '\t' + str(di[annotation][ccre][comparison]['percent']) + '\t' +
					di[annotation][ccre][comparison]['number'] + '\t' + str(di[annotation][ccre][comparison]['position']) + '\n')

