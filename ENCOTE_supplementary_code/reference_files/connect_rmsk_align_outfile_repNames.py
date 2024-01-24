'''
Connects repeat names from RepeatMasker out and align files.
Usage: python3 connect_rmsk_align_outfile_repNames.py <rmsk out file> <rmsk align file> <output file>
'''

import sys

if len(sys.argv) != 4:
	sys.exit(__doc__)

align_di = {}
with open(sys.argv[2], 'r') as f:
	for line in f:
		if line.strip() == '': #empty line, skip
			continue
		elif line.startswith(' '): #line starts with space, skip
			continue
		elif line.startswith('C'): #line starts with "C" => complementary strand of alignment, skip
			continue
		elif line.startswith('Matrix'): #line with Matrix info (whatever that is), skip
			continue
		elif line.startswith('Kimura'): #line with Kimura divergence, skip
			continue
		elif line.startswith('Transitions'): #line with transition/transversion rate/numbers, skip
			continue
		elif line.startswith('Gap_init'): #line with gap initiation rate and gap size, skip
			continue
		fields = line.rstrip('\n').split()
		score = int(fields[0])
		name = fields[4] + ':' + fields[5] + '-' + fields[6]
		if fields[8] == 'C': #See if repeat is on negative strand, get repeat name/family accordingly
			repeat_info = fields[9]
		else:
			repeat_info = fields[8]
		repeat = repeat_info.split('#')[0]
		rmsk_id = fields[-1]
		if name not in align_di:
			align_di[name] = [repeat, rmsk_id, score]
		else:
			prev_score = align_di[name][2]
			if score > prev_score: #Only replace previous annotation if the current (SW) score is greater
				align_di[name] = [repeat, rmsk_id, score]

connect_di = {}
with open(sys.argv[1], 'r') as f:
	#Skip first 3 lines of header
	f.readline()
	f.readline()
	f.readline()
	#Get rmsk output info
	for line in f:
		fields = line.rstrip('\n').split()
		name = fields[4] + ':' + fields[5] + '-' + fields[6]
		repeat = fields[9]
		rmsk_id = fields[-1]
		if name in align_di: #Only consider elements with the exact same chr,start,stop coordinates (i.e. "high confidence" and non-fragmented elements)
			align_repeat = align_di[name][0]
			align_id = align_di[name][1]
			if align_id != rmsk_id: #Check to make sure align ID is the same as out file ID
				continue #Skip if not the same for some reason
			if repeat not in connect_di:
				connect_di[repeat] = {}
				connect_di[repeat]['all'] = 0
			if align_repeat not in connect_di[repeat]:
				connect_di[repeat][align_repeat] = 0
			connect_di[repeat][align_repeat] += 1
			connect_di[repeat]['all'] += 1

with open(sys.argv[-1], 'w') as o:
	o.write('Repeat\tAlign_repeat\tNumber\tPercentage\n')
	for repeat in connect_di:
		for align_repeat in connect_di[repeat]:
			if align_repeat == 'all':
				continue
			o.write(repeat + '\t' + align_repeat + '\t' + str(connect_di[repeat][align_repeat]) + '\t' + str(connect_di[repeat][align_repeat]/connect_di[repeat]['all']) + '\n')

