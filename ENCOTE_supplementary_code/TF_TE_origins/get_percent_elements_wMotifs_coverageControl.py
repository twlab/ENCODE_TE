'''
Gets the percentage of elements in the subfamily with each motif.
This version is for coverage controlled and length controlled shared motifs.
Outputs two files:
	1. consensus controlled background percentages
	2. length controlled background percentages
Usage: python3 get_percent_elements_wMotifs_coverageControl.py <directory of subfamily confident motifs> <kimura divergence file> <consensus length file> <output file basename>
'''

import sys, os

if len(sys.argv) != 5:
	sys.exit(__doc__)

dir_files = sys.argv[1].rstrip('/') + '/'
try:
	list_files = os.scandir(path=dir_files)
except OSError:
	sys.stderr.write(sys.argv[1] + ' directory not found. Exiting...')
	raise
	sys.exit(__doc__)

list_classes = ["DNA", "LTR", "SINE", "LINE", "SVA", "ERV-int", "Other"]
TE_classes = {'subfamily': {}, 'class': {}}
kimura_div = {}
with open(sys.argv[2], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		divergence = fields[2]
		if "DNA" in fields[1]:
			TEclass = "DNA"
		elif "LTR" in fields[1]:
			if "-int" in subfamily or "MamGypsy2" in subfamily:
				TEclass = "ERV-int"
			else:
				TEclass = "LTR"
		elif "SINE" in fields[1]:
			TEclass = "SINE"
		elif "LINE" in fields[1]:
			TEclass = "LINE"
		elif "SVA" in fields[1]:
			TEclass = "SVA"
		else:
			TEclass = "Other"
		TE_classes['subfamily'][subfamily] = TEclass
		kimura_div[subfamily] = divergence

lengths = {}
with open(sys.argv[3], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		length = fields[1]
		lengths[subfamily] = length

di = {}
for filename in list_files:
	with open(filename.path, 'r') as f:
		header = f.readline()
		for line in f:
			fields = line.rstrip('\n').split('\t')
			subfamily = fields[0]
			motif = fields[1]
			cluster_name = fields[12]
			sig_motifs = fields[13]
			#numbers from coverage controlled background
			total = int(fields[5]) + int(fields[6]) + int(fields[7]) + int(fields[8])
			wMotif = int(fields[5]) + int(fields[6])
			wMotif_wAnnotation = int(fields[5])
			wMotif_noAnnotation = int(fields[6])
			noMotif_wAnnotation = int(fields[7])
			noMotif_noAnnotation = int(fields[8])
			percent = round(wMotif/total*100, 2)
			percent_wAnnotation = round(wMotif_wAnnotation/(wMotif_wAnnotation+noMotif_wAnnotation)*100, 2)
			percent_noAnnotation = round(wMotif_noAnnotation/(wMotif_noAnnotation+noMotif_noAnnotation)*100, 2)
			if subfamily not in di:
				di[subfamily] = {}
			if subfamily in TE_classes['subfamily']:
				TEclass = TE_classes['subfamily'][subfamily]
			else:
				TEclass = 'NA'
			if subfamily in kimura_div:
				divergence = kimura_div[subfamily]
			else:
				divergence = 'NA'
			if subfamily in lengths:
				length = lengths[subfamily]
			else:
				length = 'NA'
			di[subfamily][motif] = {}
			di[subfamily][motif]['coverage'] = [subfamily, TEclass, divergence, length, cluster_name, motif, str(percent), str(percent_wAnnotation), str(percent_noAnnotation), sig_motifs]
			#numbers from coverage controlled background
			total = int(fields[14]) + int(fields[15]) + int(fields[16]) + int(fields[17])
			wMotif = int(fields[14]) + int(fields[15])
			wMotif_wAnnotation = int(fields[14])
			wMotif_noAnnotation = int(fields[15])
			noMotif_wAnnotation = int(fields[16])
			noMotif_noAnnotation = int(fields[17])
			percent = round(wMotif/total*100, 2)
			percent_wAnnotation = round(wMotif_wAnnotation/(wMotif_wAnnotation+noMotif_wAnnotation)*100, 2)
			percent_noAnnotation = round(wMotif_noAnnotation/(wMotif_noAnnotation+noMotif_noAnnotation)*100, 2)
			di[subfamily][motif]['length'] = [subfamily, TEclass, divergence, length, cluster_name, motif, str(percent), str(percent_wAnnotation), str(percent_noAnnotation), sig_motifs]

with open(sys.argv[-1] + '_coverageControl', 'w') as o:
	o.write('Subfamily\tClass\tKimura_div\tLength\tMotif_cluster\tMotif\tPercent_motifs\tPercent_cCRE\tPercent_nonCRE\tSig_motifs\n')
	for subfamily in sorted(di):
		for motif in di[subfamily]:
			o.write('\t'.join(di[subfamily][motif]['coverage']) + '\n')

with open(sys.argv[-1] + '_lengthControl', 'w') as o:
	o.write('Subfamily\tClass\tKimura_div\tLength\tMotif_cluster\tMotif\tPercent_motifs\tPercent_cCRE\tPercent_nonCRE\tSig_motifs\n')
	for subfamily in sorted(di):
		for motif in di[subfamily]:
			o.write('\t'.join(di[subfamily][motif]['length']) + '\n')

