'''
Creates the TE subfamily-cCRE enrichment matrix for plotting in R
Add in class info for TE subfamilies
Usage: python3 create_subfamily_cCRE_enrichment_matrix_forR.py <TE-cCRE enrichment file> <TE-class association file> <file with TE subfamilies to include> <"count" or "bases"> <output file>
'''

import sys

if len(sys.argv) != 6:
	sys.exit(__doc__)

enrich_type = sys.argv[4]
if enrich_type != "count" and enrich_type != "bases":
	print('<"count" or "bases"> must be one of the two options')
	sys.exit(__doc__)

enrich_di = {}
enrich_ccres = []
with open(sys.argv[1], 'r') as f:
	header = f.readline().strip().split('\t')
	select_indices = []
	for i in range(len(header)):
		if enrich_type in header[i]:
			select_indices.append(i)
			enrich_ccres.append(header[i].split('(')[0])
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		enrich_scores = []
		for index in select_indices:
			enrich_scores.append(fields[index])
		enrich_di[subfamily] = enrich_scores

class_di = {}
with open(sys.argv[2], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		if "DNA" in fields[1]:
			te_class = "DNA"
		elif "LTR" in fields[1]:
			te_class = "LTR"
		elif "LINE" in fields[1]:
			te_class = "LINE"
		elif "SINE" in fields[1]:
			te_class = "SINE"
		elif "SVA" in fields[1]:
			te_class = "SVA"
		else:
			te_class = "Other"
		class_di[subfamily] = te_class

included_subfamilies = []
with open(sys.argv[3], 'r') as f:
	for line in f:
		included_subfamilies.append(line.strip())

te_classes = ['DNA', 'LTR', 'SINE', 'LINE', 'SVA', 'Other']
with open(sys.argv[-1], 'w') as o:
	o.write('Class\tSubfamily\tcCRE\tFoldEnrich\n')
	#Looks weird, but this is to group the subfamilies by TE class
	for te_class in te_classes:
		for subfamily in included_subfamilies:
			if subfamily in class_di:
				if class_di[subfamily] == te_class:
					for i in range(len(enrich_ccres)):
						o.write(te_class + '\t' + subfamily + '\t' + enrich_ccres[i] + '\t' + enrich_di[subfamily][i] + '\n')

