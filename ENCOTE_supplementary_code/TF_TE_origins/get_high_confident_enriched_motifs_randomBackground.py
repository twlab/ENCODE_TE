'''
Gets high confident enriched motifs based on Fisher's exact test results.
Outputs median of Fisher's exact tests (based on odds ratio).
Usage: python3 get_high_confident_enriched_motifs_randomBackground.py <enriched motifs AME info file> <enriched motifs fishers directory> <output file>
'''

import sys, os
import statistics as stats

#Fisher's exact test thresholds for high confidence enriched motifs
p_thresh = 0.05 #Fisher's exact test p-value
odds_thresh = 2 #Odds ratio
true_pos_thresh = 10 #Number of +Motif+Annotation elements

if len(sys.argv) != 4:
	sys.exit(__doc__)

ame_info = {}
ame_motifs = []
with open(sys.argv[1], 'r') as f:
	header = f.readline().rstrip('\n').split('\t')
	ame_header = '\t'.join(header[:2]) + '\t' + header[-1]
	for line in f:
		fields = line.rstrip('\n').split('\t')
		motif = fields[2]
		info = fields[:2]
		info.append(fields[-1])
		ame_info[motif] = info
		ame_motifs.append(motif)

filenames = []
output_dir = sys.argv[2].rstrip('/') + '/'
for filename in os.scandir(output_dir):
	if filename.is_file():
		filenames.append(filename.path)

fishers_info = {}
passing_motifs = {}
for filename in filenames:
	with open(filename, 'r') as f:
		fishers_header = f.readline().rstrip('\n').split('\t')
		fishers_header[3] = 'FishersP'
		for line in f:
			fields = line.rstrip('\n').split('\t')
			motif = fields[1]
			info = fields
			odds = float(fields[2])
			p_val = float(fields[3])
			true_pos = int(fields[4])
			perc_motif_annotation = round(int(fields[4])/(int(fields[4]) + int(fields[6]))*100)
			perc_motif_nonAnnotation = round(int(fields[5])/(int(fields[5]) + int(fields[7]))*100)
			info.append(str(perc_motif_annotation))
			info.append(str(perc_motif_nonAnnotation))
			if motif not in fishers_info:
				fishers_info[motif] = {}
				fishers_info[motif]['odds'] = [[], []]
			fishers_info[motif][filename] = info
			fishers_info[motif]['odds'][0].append(odds)
			fishers_info[motif]['odds'][1].append(filename)
			if motif in ame_motifs:
				if odds >= odds_thresh:
					if p_val <= p_thresh:
						if true_pos >= true_pos_thresh:
							if motif not in passing_motifs:
								passing_motifs[motif] = 0
							passing_motifs[motif] += 1

with open(sys.argv[-1], 'w') as o:
	o.write('\t'.join(fishers_header) + '\tPerc_Annotation+Motif\tPerc_Annotation-Motif\t' + ame_header + '\n')
	for motif in ame_motifs:
		if motif in passing_motifs:
			if passing_motifs[motif] == len(filenames): #All randomizations have motif as enriched by Fisher's exact test criteria
				median_odds = stats.median_low(fishers_info[motif]['odds'][0])
				name_index = fishers_info[motif]['odds'][0].index(median_odds)
				filename = fishers_info[motif]['odds'][1][name_index]
				o.write('\t'.join(fishers_info[motif][filename]) + '\t' + '\t'.join(ame_info[motif]) + '\n')
			else:
				continue

