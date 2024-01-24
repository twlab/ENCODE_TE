'''
Calculate TE proportions according to our re-annotated TE file
Usage: python3 calc_TE_proportions.py <input reAnno bed file> <output file>
'''

import sys

genome_size = 3209286105 #~3.2 billion bp genome size, number from adding up chromosome sizes from /bar/genomes/hg38/hg38.chrom.sortV.sizes

if len(sys.argv) != 3:
	sys.exit(__doc__)

di = {}
with open(sys.argv[1], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		TE_len = int(fields[2]) - int(fields[1])
		if "DNA" in fields[4]:
			TEclass = "DNA"
		elif "LINE" in fields[4]:
			TEclass = "LINE"
		elif "SINE" in fields[4]:
			TEclass = "SINE"
		elif "LTR" in fields[4]:
			TEclass = "LTR"
		elif "SVA" in fields[4]:
			TEclass = "SVA"
		else:
			TEclass = "Other"
		if TEclass not in di:
			di[TEclass] = 0
		di[TEclass] += TE_len

with open(sys.argv[-1], 'w') as o:
	o.write('Class\tGenomic_proportion\tGenomic_percentage\n')
	allTE = 0
	for TEclass in di:
		if TEclass != "Other":
			o.write(TEclass + '\t' + str(di[TEclass]/genome_size) + '\t' + str(round(di[TEclass]/genome_size*100, 2)) + '%\n')
			allTE += di[TEclass]
	o.write('non-TE\t' + str((genome_size - allTE)/genome_size) + '\t' + str(round((genome_size - allTE)/genome_size*100, 2)) + '%\n')

