'''
Gets high confident enriched motifs based on Fisher's exact test results.
Combines AME and Fisher's exact test info.
Usage: python3 get_high_confident_enriched_motifs.py <enriched motifs AME info file> <enriched motifs fishers file> <output file>
'''

import sys

#Fisher's exact test thresholds for high confidence enriched motifs
p_thresh = 0.05 #Fisher's exact test p-value
odds_thresh = 2 #Odds ratio
true_pos_thresh = 10 #Number of +Motif+Annotation elements

if len(sys.argv) != 4:
	sys.exit(__doc__)

ame_info = {}
ame_motifs = []
with open(sys.argv[1], 'r') as f:
	ame_header = f.readline().rstrip('\n').split('\t')
	for line in f:
		fields = line.rstrip('\n').split('\t')
		motif = fields[2]
		info = fields[:2] + fields[3:]
		ame_info[motif] = info
		ame_motifs.append(motif)

fishers_info = {}
passing_motifs = []
with open(sys.argv[2], 'r') as f:
	fishers_header = f.readline().rstrip('\n').split('\t')
	fishers_header[3] = 'FishersP'
	for line in f:
		fields = line.rstrip('\n').split('\t')
		motif = fields[1]
		info = fields
		perc_motif_annotation = int(fields[4])/(int(fields[4]) + int(fields[6]))*100
		perc_motif_nonAnnotation = int(fields[5])/(int(fields[5]) + int(fields[7]))*100
		info.append(str(perc_motif_annotation))
		info.append(str(perc_motif_nonAnnotation))
		fishers_info[motif] = info
		odds = float(fields[2])
		p_val = float(fields[3])
		true_pos = int(fields[4])
		if motif in ame_motifs:
			if odds >= odds_thresh:
				if p_val <= p_thresh:
					if true_pos >= true_pos_thresh:
						passing_motifs.append(motif)

with open(sys.argv[-1], 'w') as o:
	o.write('\t'.join(fishers_header) + '\tPerc_Annotation+Motif\tPerc_Annotation-Motif\t' + '\t'.join(ame_header[:2] + ame_header[3:]) + '\n')
	for motif in passing_motifs:
		o.write('\t'.join(fishers_info[motif]) + '\t' + '\t'.join(ame_info[motif]) + '\n')

