'''
Combines log2 enrichments for TE subfamilies in each cCRE type across full classified cell/tissue types.
Outputs each cCRE into separate files
Takes input for the minimum number of overlaps required to report log2 enrichment
TE subfamily is only reported if at least one enrichment passes the overlap threshold
Usage: python3 combine_enrichments_TEsubfamilies_fullClassCCREs_forR.py <TE subfamily enrichment directory> <TE subfamily-class association file> <min number of overlaps required> <output file basename>
'''

import sys, os

if len(sys.argv) != 5:
	sys.exit(__doc__)

try:
	min_overlaps = int(sys.argv[3])
except TypeError:
	print('<min number of overlaps required> must be an integer')
	sys.exit(__doc__)

dir_path = sys.argv[1].rstrip('/') + '/'
list_files = os.listdir(dir_path)
ccre_types = []
subfamilies = []
celltypes = []
di = {}
include_subfamilies = []
for filename in list_files:
	with open(dir_path + filename, 'r') as f:
		header = f.readline()
		for line in f:
			fields = line.rstrip('\n').split('\t')
			subfamily = fields[0]
			if subfamily not in subfamilies:
				subfamilies.append(subfamily)
			celltype = fields[1]
			if celltype not in celltypes:
				celltypes.append(celltype)
			ccre = fields[2]
			if ccre not in ccre_types:
				ccre_types.append(ccre)
			num_overlaps = int(fields[4])
			enrichment = fields[5]
			key = subfamily + '/' + celltype + '/' + ccre
			if num_overlaps >= min_overlaps:
				di[key] = enrichment
				if subfamily not in include_subfamilies:
					include_subfamilies.append(subfamily)
			else:
				di[key] = enrichment

list_classes = ["DNA", "LTR", "SINE", "LINE", "SVA"]
di_class = {}
with open(sys.argv[2], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		if "DNA" in fields[1]:
			TEclass = "DNA"
		elif "LTR" in fields[1]:
			TEclass = "LTR"
		elif "SINE" in fields[1]:
			TEclass = "SINE"
		elif "LINE" in fields[1]:
			TEclass = "LINE"
		elif "SVA" in fields[1]:
			TEclass = "SVA"
		else:
			continue
		di_class[subfamily] = TEclass

for ccre in ccre_types:
	filename = sys.argv[-1] + '_' + ccre
	with open(filename, 'w') as o:
		o.write('Subfamily\tClass\tCelltype\tcCRE\tEnrichment\n')
		for subfamily in include_subfamilies:
			for celltype in celltypes:
				key = subfamily + '/' + celltype + '/' + ccre
				if subfamily in di_class:
					o.write(subfamily + '\t' + di_class[subfamily] + '\t' + celltype + '\t' + ccre + '\t' + di[key] + '\n')

