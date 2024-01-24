'''
Gets enriched motif annotations for each TE of a subfamily from FIMO motif scan.
Assumes all elements are from the same subfamily.
Usage: python3 get_enriched_motif_annotations.py <input subfamily bed file> <path to dir of TE subfamily FIMO outputs> <output file>
'''

import sys, os

if len(sys.argv) != 4:
	sys.exit(__doc__)

list_elements = []
with open(sys.argv[1], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[3]
		name = fields[0] + ':' + fields[1] + '-' + fields[2] + '(' + fields[5] + ')'
		list_elements.append(name)

dir_path = sys.argv[2].rstrip('/') + '/'
motif_scans = os.listdir(dir_path)

di_motifs = {}
motif_list = []
for motif_scan in motif_scans:
	filename = dir_path + motif_scan + '/fimo.tsv'
	with open(filename, 'r') as f:
		header = f.readline()
		for line in f:
			if line.startswith('#'):
				continue
			if line.strip() == '':
				continue
			fields = line.rstrip('\n').split('\t')
			if fields[0] == '':
				motif = fields[1]
			else:
				motif = fields[0]
			if motif not in motif_list:
				motif_list.append(motif)
			name = fields[2]
			if name not in di_motifs:
				di_motifs[name] = []
			while len(di_motifs[name]) < len(motif_list):
				di_motifs[name].append(0)
			motif_index = motif_list.index(motif)
			di_motifs[name][motif_index] += 1

for name in list_elements:
	if name not in di_motifs:
		di_motifs[name] = []
	while len(di_motifs[name]) < len(motif_list):
		di_motifs[name].append(0)

with open(sys.argv[-1], 'w') as o:
	#Header
	o.write('Subfamily\tElement')
	for motif in motif_list:
		o.write('\t' + motif)
	o.write('\n')
	#Write enriched motif annotations for each TE
	for name in list_elements:
		if len(di_motifs[name]) == 0:
			o.write(subfamily + '\t' + name + '\n')
		else:
			o.write(subfamily + '\t' + name + '\t' + '\t'.join([str(x) for x in di_motifs[name]]) + '\n')

