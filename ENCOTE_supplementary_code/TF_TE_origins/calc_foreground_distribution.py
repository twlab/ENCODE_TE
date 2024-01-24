'''
Calculates distribution of TEs for foreground/background elements relative to consensus.
Outputs all distribution coverages for foreground/background to <output file basename> separated by TE subfamily.
Decides which motifs are likely to be enriched simply by virtue of consensus distribution. Prints these motifs to standard output.
Usage: python3 calc_foreground_distribution.py <list of TE subfamilies> <TE consensus motifs file> <directory of overlap/non-overlap files> <directory of alignFastas> <output file basename>
'''

import sys, os

exclude_thresh = 2 #Threshold for minimum foreground/background enrichment in coverage required to exclude motifs in the region
exclude_base = 10 #Threshold for minimum foreground coverage required to exclude motifs in the region

if len(sys.argv) != 6:
	sys.exit(__doc__)

list_subfamilies = []
with open(sys.argv[1], 'r') as f:
	for line in f:
		subfamily = line.rstrip('\n')
		list_subfamilies.append(subfamily)

consensus_motifs = {}
with open(sys.argv[2], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		if subfamily not in consensus_motifs:
			consensus_motifs[subfamily] = {}
		motif = fields[1]
		starts = fields[2].split(',')
		stops = fields[3].split(',')
		consensus_motifs[subfamily][motif] = []
		for i in range(len(starts)):
			start = int(starts[i]) - 1 #Convert to 0-based
			stop = int(stops[i]) - 1 #Convert to 0-based
			consensus_motifs[subfamily][motif].append([start, stop])

overlap_dir = sys.argv[3].rstrip('/') + '/'
all_overlap_files = os.listdir(overlap_dir)
di_overlap = {}
di_nonOverlap = {}
for subfamily in list_subfamilies:
	overlap_file = subfamily + '.overlap'
	nonOverlap_file = subfamily + '.nonOverlap'
	if overlap_file in all_overlap_files:
		di_overlap[subfamily] = {}
		di_overlap[subfamily]['total'] = 0
		with open(overlap_dir + overlap_file, 'r') as f:
			for line in f:
				fields = line.rstrip('\n').split('\t')
				name = fields[0] + ':' + fields[1] + '-' + fields[2] + '(' + fields[5] + ')'
				di_overlap[subfamily][name] = 'Y'
				di_overlap[subfamily]['total'] += 1
	if nonOverlap_file in all_overlap_files:
		di_nonOverlap[subfamily] = {}
		di_nonOverlap[subfamily]['total'] = 0
		with open(overlap_dir + nonOverlap_file, 'r') as f:
			for line in f:
				fields = line.rstrip('\n').split('\t')
				name = fields[0] + ':' + fields[1] + '-' + fields[2] + '(' + fields[5] + ')'
				di_nonOverlap[subfamily][name] = 'N'
				di_nonOverlap[subfamily]['total'] += 1

align_dir = sys.argv[4].rstrip('/') + '/'
list_align_files = os.listdir(align_dir)
coverage_di = {}
lengths_di = {}
for subfamily in list_subfamilies:
	filename = subfamily + '.alignFasta'
	if filename in list_align_files:
		if subfamily in di_overlap:
			coverage_di[subfamily] = {}
			lengths_di[subfamily] = {}
			with open(align_dir + filename, 'r') as f:
				for line in f:
					if line.startswith('>'):
						name = line.lstrip('>').rstrip('\n')
						seq = f.readline().rstrip('\n')
						ref = f.readline().rstrip('\n')
						coverage = []
						for i in range(len(seq)):
							if seq[i] == '-': #Deletion in element => no coverage at reference base
								coverage.append(0)
							else:
								if ref[i] == '-': #Insertion in element => no reference base
									continue
								else: #Base in element and reference
									coverage.append(1)
						if name in di_overlap[subfamily]:
							if 'overlap' not in coverage_di[subfamily]:
								coverage_di[subfamily]['overlap'] = coverage
								lengths_di[subfamily]['overlap'] = {name: str(len(seq.replace('-', '')))}
							else:
								for i in range(len(coverage)):
									coverage_di[subfamily]['overlap'][i] += coverage[i]
								lengths_di[subfamily]['overlap'][name] = str(len(seq.replace('-', '')))
						elif name in di_nonOverlap[subfamily]:
							if 'nonOverlap' not in coverage_di[subfamily]:
								coverage_di[subfamily]['nonOverlap'] = coverage
								lengths_di[subfamily]['nonOverlap'] = {name: str(len(seq.replace('-', '')))}
							else:
								for i in range(len(coverage)):
									coverage_di[subfamily]['nonOverlap'][i] += coverage[i]
								lengths_di[subfamily]['nonOverlap'][name] = str(len(seq.replace('-', '')))
						else: #Should not get here
							sys.stderr.write('something is wrong\n')
							sys.exit()

di_exclude_motifs = {}
for subfamily in list_subfamilies:
	if subfamily in coverage_di:
		if 'overlap' in coverage_di[subfamily] and 'nonOverlap' in coverage_di[subfamily]:
			exclude_coverage = []
			enrich_coverage = []
			for i in range(len(coverage_di[subfamily]['overlap'])):
				try: #Enrichment calculation = relative coverage in foreground (elements with overlap divided by total number of elements) divided by relative coverage in background
					enrich = (coverage_di[subfamily]['overlap'][i] / di_overlap[subfamily]['total']) / (coverage_di[subfamily]['nonOverlap'][i] / di_nonOverlap[subfamily]['total'])
				except ZeroDivisionError:
					if di_overlap[subfamily]['total'] == 0: #If error caused by 0 in denominator of foreground (no elements in foreground), set enrichment to 0
						enrich = 0
					elif coverage_di[subfamily]['overlap'][i] == 0 and coverage_di[subfamily]['nonOverlap'][i] == 0: #Both foreground and background have no coverage
						enrich = 1
					else: #If error caused by 0 in background (no coverage or elements in background), set enrichment to max
						enrich = 100
				if coverage_di[subfamily]['overlap'][i] >= exclude_base and enrich >= exclude_thresh:
					exclude_coverage.append('Y')
				else:
					exclude_coverage.append('N')
				enrich_coverage.append(enrich)
			for motif in consensus_motifs[subfamily]:
				for start, stop in consensus_motifs[subfamily][motif]:
					regional_coverage = exclude_coverage[start:(stop+1)]
					if 'N' in regional_coverage: #If any parts of the motif region are not different between foreground and background, skip motif
						continue
					else: #All consensus base positions of motif region are found more often in foreground than background, add motif to exclusion list
						if subfamily not in di_exclude_motifs:
							di_exclude_motifs[subfamily] = []
						di_exclude_motifs[subfamily].append([motif, str(start), str(stop)])
			with open(sys.argv[-1] + '_' + subfamily, 'w') as o:
				o.write('foreground\t' + '\t'.join([str(x) for x in coverage_di[subfamily]['overlap']]) + '\n')
				o.write('background\t' + '\t'.join([str(x) for x in coverage_di[subfamily]['nonOverlap']]) + '\n')
				o.write('enrichment\t' + '\t'.join([str(x) for x in enrich_coverage]) + '\n')

for subfamily in list_subfamilies:
	if subfamily in lengths_di:
		with open(sys.argv[-1] + '_' + subfamily + '.lengths', 'w') as o:
			o.write('Element\tAnnotation\tLength\n')
			if 'overlap' in lengths_di[subfamily]:
				for name in lengths_di[subfamily]['overlap']:
					o.write(name + '\tOverlap\t' + lengths_di[subfamily]['overlap'][name] + '\n')
			if 'nonOverlap' in lengths_di[subfamily]:
				for name in lengths_di[subfamily]['nonOverlap']:
					o.write(name + '\tNon-overlap\t' + lengths_di[subfamily]['nonOverlap'][name] + '\n')

print('Subfamily\tMotif\tStart\tStop')
for subfamily in list_subfamilies:
	if subfamily in di_exclude_motifs:
		for i in range(len(di_exclude_motifs[subfamily])):
			print(subfamily + '\t' + '\t'.join(di_exclude_motifs[subfamily][i]))

