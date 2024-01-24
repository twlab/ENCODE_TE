'''
Creates the TE class-cCRE enrichment matrix for plotting in R
Usage: python3 create_class_cCRE_enrichment_matrix_forR.py <TE class-cCRE enrichment file> <"count" or "bases"> <output file>
'''

import sys

if len(sys.argv) != 4:
	sys.exit(__doc__)

enrich_type = sys.argv[2]
if enrich_type != "count" and enrich_type != "bases":
	print('<"count" or "bases"> must be one of the two options')
	sys.exit(__doc__)

enrich_di = {}
enrich_ccres = []
te_classes = []
with open(sys.argv[1], 'r') as f:
	header = f.readline().strip().split('\t')
	select_indices = []
	for i in range(len(header)):
		if enrich_type in header[i]:
			select_indices.append(i)
			enrich_ccres.append(header[i].split('(')[0])
	for line in f:
		fields = line.rstrip('\n').split('\t')
		te_class = fields[0]
		enrich_scores = []
		for index in select_indices:
			enrich_scores.append(fields[index])
		enrich_di[te_class] = enrich_scores
		te_classes.append(te_class)

with open(sys.argv[-1], 'w') as o:
	o.write('Class\tcCRE\tFoldEnrich\n')
	for te_class in te_classes:
		for i in range(len(enrich_ccres)):
			o.write(te_class + '\t' + enrich_ccres[i] + '\t' + enrich_di[te_class][i] + '\n')

