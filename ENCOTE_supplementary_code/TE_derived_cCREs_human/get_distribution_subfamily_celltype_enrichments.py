'''
Gets the distribution of TE subfamilies that are enriched for each cCRE for each TE class
Usage: python3 get_distribution_subfamily_celltype_enrichments.py <input forR enrichment file with class info> <log2 enrichment threshold> <output file>
'''

import sys

if len(sys.argv) != 4:
	sys.exit(__doc__)

try:
	thresh = float(sys.argv[2])
except TypeError:
	print('<log2 enrichment threshold> must be a real number')
	sys.exit(__doc__)

celltypes = []
TEclasses = []
ccre_types = []
di = {}
di_assoc = {}
with open(sys.argv[1], 'r') as f:
	for line in f:
		if '=>' in line:
			continue
		elif line.rstrip() == '':
			continue
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		TEclass = fields[1]
		if TEclass not in TEclasses:
			TEclasses.append(TEclass)
		if subfamily not in di_assoc:
			di_assoc[subfamily] = TEclass
		celltype = fields[2]
		if celltype not in celltypes:
			celltypes.append(celltype)
		ccre = fields[3]
		if ccre not in ccre_types:
			ccre_types.append(ccre)
		key = subfamily + '/' + ccre
		if key not in di:
			di[key] = 0
		if fields[4] != 'NA':
			enrich = float(fields[4])
			if enrich >= thresh:
				di[key] += 1
		else:
			continue

with open(sys.argv[-1], 'w') as o:
	o.write('Class\tcCRE\tNumEnrichedCelltypes\tNumSubfamilies\n')
	for TEclass in TEclasses:
		for ccre in ccre_types:
			numbers = []
			for i in range(len(celltypes)+1):
				numbers.append(0)
			for subfamily in di_assoc:
				if di_assoc[subfamily] == TEclass:
					key = subfamily + '/' + ccre
					index = di[key]
					numbers[index] += 1
			for i in range(len(numbers)):
				o.write(TEclass + '\t' + ccre + '\t' + str(i) + '\t' + str(numbers[i]) + '\n')

