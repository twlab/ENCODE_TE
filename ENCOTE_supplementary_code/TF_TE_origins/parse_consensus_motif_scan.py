'''
Parse FIMO output for motif scan of TE subfamily consensus sequences.
Assumes motif is in the first column of the FIMO output (like when run with HOCOMOCO database).
Usage: python3 parse_consensus_motif_scan.py <input fimo tsv> <list of TE subfamilies file> <output file>
'''

import sys

if len(sys.argv) != 4:
	sys.exit(__doc__)

subfamilies = []
with open(sys.argv[2], 'r') as f:
	for line in f:
		subfamily = line.strip()
		subfamilies.append(subfamily)

di = {}
motifs = []
with open(sys.argv[1], 'r') as f:
	header = f.readline()
	for line in f:
		if line.startswith('#'):
			continue
		elif line.strip() == '':
			continue
		fields = line.rstrip('\n').split('\t')
		motif = fields[0]
		subfamily = fields[2]
		info = fields[3:6]
		if motif not in motifs:
			motifs.append(motif)
		index = motifs.index(motif)
		if subfamily not in di:
			di[subfamily] = []
		while len(di[subfamily]) < len(motifs):
			di[subfamily].append([])
		if di[subfamily][index] == []: #no same motif found yet
			di[subfamily][index].append(motif)
			di[subfamily][index] += info
		else:
			for i in range(len(info)):
				di[subfamily][index][i+1] += ',' + info[i]

for subfamily in subfamilies:
	if subfamily not in di:
		continue
	while len(di[subfamily]) < len(motifs):
		di[subfamily].append([])

with open(sys.argv[-1], 'w') as o:
	o.write('Subfamily\tMotif\tStart\tStop\tStrand\n')
	for subfamily in subfamilies:
		if subfamily not in di:
			continue
		for i in range(len(motifs)):
			motif = motifs[i]
			if di[subfamily][i] != []:
				o.write(subfamily + '\t' + '\t'.join(di[subfamily][i]) + '\n')

