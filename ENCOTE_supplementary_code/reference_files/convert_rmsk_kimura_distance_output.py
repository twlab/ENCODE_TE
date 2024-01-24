'''
Converts Kimura distances from RepeatMasker "calcDivergenceFromAlign.pl" to a more python friendly format.
Tries the following to get the best Kimura % divergence for TE subfamily names that do not match rmsk out names:
	1. If there is an "-int" in the rmsk out name, check if removing the "-int" creates a recognizable subfamily name from the align file.
	2. If the json file has 3' end annotation or a single annotation for the subfamily, use the 3' or single annotation.
	3. Else, use the most abundant (by copy) annotation
Prints TE subfamily-rmsk out name connections to standard output.
Usage: python3 convert_rmsk_kimura_distance_output.py <input file> <list of TE subfamilies file> <hg38 LINE json file> <output file>
'''

import sys
import json

if len(sys.argv) != 5:
	sys.exit(__doc__)

subfamilies = []
with open(sys.argv[2], 'r') as f:
	for line in f:
		subfamilies.append(line.strip())

kimura_input = {}
with open(sys.argv[1], 'r') as f:
	for i in range(4): #Skip the first 4 lines
		f.readline()
	#Get header
	header = f.readline().split()
	new_header = header[1] + '\t' + header[0] + '\t' + header[-1] + '\n'
	#Skip next line
	f.readline()
	#Get all of the Kimura % divergences
	for line in f:
		if line.strip() == '': #Reached the end of the Kimura % divergences
			break
		fields = line.rstrip().split()
		kimura_input[fields[1]] = [fields[0], fields[-1]]

connect_di = {}
with open(sys.argv[3], 'r') as json_file:
	json_info = json.load(json_file)
	for subfamily in json_info:
		if '3end' in json_info[subfamily]:
			copies = []
			positions = []
			for position in list(json_info[subfamily]):
				positions.append(position)
				copies.append(int(json_info[subfamily][position]['cp']))
			if int(json_info[subfamily]['3end']['cp']) >= (max(copies)*0.2): #If the 3'end has at least 20% of the max number of copies, use it for kimura age
				name = json_info[subfamily]['3end']['te']
				connect_di[subfamily] = name
			else: #identified 3'end is unlikely to be real 3'end, use the most abundant piece
				index = copies.index(max(copies))
				position = positions[index]
				name = json_info[subfamily][position]['te']
				connect_di[subfamily] = name
		elif len(json_info[subfamily]) == 1:
			position = list(json_info[subfamily])[0]
			name = json_info[subfamily][position]['te']
			connect_di[subfamily] = name
		else:
			copies = []
			positions = []
			for position in list(json_info[subfamily]):
				positions.append(position)
				copies.append(int(json_info[subfamily][position]['cp']))
			index = copies.index(max(copies))
			position = positions[index]
			name = json_info[subfamily][position]['te']
			connect_di[subfamily] = name

with open(sys.argv[-1], 'w') as o:
	o.write(new_header)
	for subfamily in subfamilies:
		if subfamily in kimura_input:
			print(subfamily + '\t' + subfamily)
			o.write(subfamily + '\t' + '\t'.join(kimura_input[subfamily]) + '\n')
		else:
			if '-int' in subfamily:
				basename = '-'.join(subfamily.split('-')[:-1])
				if basename in kimura_input:
					o.write(subfamily + '\t' + '\t'.join(kimura_input[basename]) + '\n')
					print(subfamily + '\t' + basename)
					continue
			if subfamily in connect_di:
				frag_name = connect_di[subfamily]
				if frag_name in kimura_input:
					o.write(subfamily + '\t' + '\t'.join(kimura_input[frag_name]) + '\n')
				else:
					print(frag_name + ' for ' + subfamily + ' does not have Kimura % divergence', file=sys.stderr)
				print(subfamily + '\t' + frag_name)
			else:
				print('Could not find Kimura % divergence match for ' + subfamily, file=sys.stderr)
				print(subfamily + '\tNA')

