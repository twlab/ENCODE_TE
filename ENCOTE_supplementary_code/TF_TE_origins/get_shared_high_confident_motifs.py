'''
Gets high confident enriched motifs based on Fisher's exact test results.
Combines high confident motifs from length controlled (AME and Fisher's exact tests) and coverage controlled (Fisher's exact test).
Usage: python3 get_shared_high_confident_motifs.py <length controlled high confident motifs file> <coverage controlled fishers file> <output file>
'''

import sys

#Fisher's exact test thresholds for high confidence enriched motifs
p_thresh = 0.05 #Fisher's exact test p-value
odds_thresh = 2 #Odds ratio
true_pos_thresh = 10 #Number of +Motif+Annotation elements

if len(sys.argv) != 4:
	sys.exit(__doc__)

length_info = {}
length_motifs = []
with open(sys.argv[1], 'r') as f:
	length_header = f.readline().rstrip('\n').split('\t')
	for line in f:
		fields = line.rstrip('\n').split('\t')
		motif = fields[1]
		length_info[motif] = fields[10:13] + fields[4:8]
		length_motifs.append(motif)

coverage_info = {}
passing_motifs = []
with open(sys.argv[2], 'r') as f:
	fishers_header = f.readline().rstrip('\n')
	for line in f:
		fields = line.rstrip('\n').split('\t')
		motif = fields[1]
		info = fields
		try:
			perc_motif_annotation = round(int(fields[5])/(int(fields[5]) + int(fields[7]))*100, 2)
			perc_motif_nonAnnotation = round(int(fields[6])/(int(fields[6]) + int(fields[8]))*100, 2)
		except ZeroDivisionError:
			continue
		info.append(str(perc_motif_annotation))
		info.append(str(perc_motif_nonAnnotation))
		coverage_info[motif] = info
		odds = float(fields[2])
		p_val = float(fields[3])
		true_pos = int(fields[5])
		if motif in length_motifs:
			if odds >= odds_thresh:
				if p_val <= p_thresh:
					if true_pos >= true_pos_thresh:
						passing_motifs.append(motif)

with open(sys.argv[-1], 'w') as o:
	o.write(fishers_header + '\tPerc_Annotation+Motif\tPerc_Annotation-Motif\t' + '\t'.join(length_header[10:13] + length_header[4:8]) + '\n')
	for motif in passing_motifs:
		o.write('\t'.join(coverage_info[motif]) + '\t' + '\t'.join(length_info[motif]) + '\n')

