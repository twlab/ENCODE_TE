'''
Calculate the cCRE proportions (overlapping TEs and not) for a specific cell type.
Usage: python3 calc_cCRE_proportions_celltype.py <input overlap file> <output file>
'''

import sys

if len(sys.argv) != 3:
	sys.exit(__doc__)

all_ccres = []
overlapping_ccres = {}
num_ccres = {}
bp_ccres = {}
num_overlap = {}
bp_overlap = {}
TEclasses = []
with open(sys.argv[1], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		ccre_id = fields[3]
		overlap = int(fields[-1])
		ccre_len = int(fields[2]) - int(fields[1])
		ccre_annotation = fields[9]
		if ccre_id not in overlapping_ccres: #If cCRE not previously seen
			if 'All cCREs' not in num_ccres:
				num_ccres['All cCREs'] = 0
				bp_ccres['All cCREs'] = 0
			num_ccres['All cCREs'] += 1
			bp_ccres['All cCREs'] += ccre_len
			ccre_groups = []
			if "dELS" in ccre_annotation:
				ccre_groups.append('dELS')
			if "pELS" in ccre_annotation:
				ccre_groups.append('pELS')
			if "PLS" in ccre_annotation:
				ccre_groups.append('PLS')
			if "DNase-H3K4me3" in ccre_annotation:
				ccre_groups.append('DNase-H3K4me3')
			if "DNase-only" in ccre_annotation:
				ccre_groups.append('DNase-only')
			if "CTCF" in ccre_annotation:
				ccre_groups.append('CTCF-bound')
			if "CTCF-only" in ccre_annotation:
				ccre_groups.append('CTCF-only')
			for ccre_group in ccre_groups:
				if ccre_group not in num_ccres:
					num_ccres[ccre_group] = 0
					bp_ccres[ccre_group] = 0
				num_ccres[ccre_group] += 1
				bp_ccres[ccre_group] += ccre_len
		if overlap == 0: #No overlap
			continue
		else:
			if "DNA" in fields[15]:
				TEclass = "DNA"
			elif "LINE" in fields[15]:
				TEclass = "LINE"
			elif "SINE" in fields[15]:
				TEclass = "SINE"
			elif "LTR" in fields[15]:
				TEclass = "LTR"
			elif "SVA" in fields[15]:
				TEclass = "SVA"
			else:
				TEclass = "Other"
			if TEclass not in TEclasses:
				TEclasses.append(TEclass)
				num_overlap[TEclass] = {}
				bp_overlap[TEclass] = {}
			if ccre_id in overlapping_ccres:
				#Check if other TE overlapping ccre is better (in terms of bp overlapping)
				prev_overlap = overlapping_ccres[ccre_id][1]
				if prev_overlap < overlap: #Current overlap is better, replace previous
					prev_class = overlapping_ccres[ccre_id][0]
					num_overlap[prev_class]['All cCREs'] -= 1
					bp_overlap[prev_class]['All cCREs'] -= prev_overlap
					num_overlap[TEclass]['All cCREs'] += 1
					bp_overlap[TEclass]['All cCREs'] += overlap
					prev_ccre_groups = overlapping_ccres[ccre_id][2]
					for ccre_group in prev_ccre_groups:
						num_overlap[prev_class][ccre_group] -= 1
						bp_overlap[prev_class][ccre_group] -= prev_overlap
						if ccre_group not in num_overlap[TEclass]:
							num_overlap[TEclass][ccre_group] = 0
							bp_overlap[TEclass][ccre_group] = 0
						num_overlap[TEclass][ccre_group] += 1
						bp_overlap[TEclass][ccre_group] += overlap
					overlapping_ccres[ccre_id] = [TEclass, overlap, prev_ccre_groups]
				else: #Current overlap is worse or same, skip
					continue
			else:
				if 'All cCREs' not in num_overlap[TEclass]:
					num_overlap[TEclass]['All cCREs'] = 0
					bp_overlap[TEclass]['All cCREs'] = 0
				num_overlap[TEclass]['All cCREs'] += 1
				bp_overlap[TEclass]['All cCREs'] += overlap
				for ccre_group in ccre_groups:
					if ccre_group not in num_overlap[TEclass]:
						num_overlap[TEclass][ccre_group] = 0
						bp_overlap[TEclass][ccre_group] = 0
					num_overlap[TEclass][ccre_group] += 1
					bp_overlap[TEclass][ccre_group] += overlap
				overlapping_ccres[ccre_id] = [TEclass, overlap, ccre_groups]

with open(sys.argv[-1], 'w') as o:
	o.write('#Total cCRE num = ' + str(num_ccres) + '\n')
	o.write('#Total cCRE bp = ' + str(bp_ccres) + '\n')
	o.write('Class\tcCRE_type\tcCRE_num\tcCRE_bp\tcCRE_num_perc\tcCRE_bp_perc\n')
	for TEclass in TEclasses:
		for ccre_type in num_overlap[TEclass]:
			o.write(TEclass + '\t' + ccre_type + '\t' + str(num_overlap[TEclass][ccre_type]) + '\t' + str(bp_overlap[TEclass][ccre_type]) + '\t' + 
				str(num_overlap[TEclass][ccre_type]/num_ccres[ccre_type]) + '\t' + str(bp_overlap[TEclass][ccre_type]/bp_ccres[ccre_type]) + '\n')
	TE_nums = {}
	TE_bps = {}
	all_nums = {}
	all_bps = {}
	for TEclass in TEclasses:
		for ccre_type in num_overlap[TEclass]:
			if ccre_type not in TE_nums:
				TE_nums[ccre_type] = 0
				TE_bps[ccre_type] = 0
			if TEclass != "Other":
				TE_nums[ccre_type] += num_overlap[TEclass][ccre_type]
				TE_bps[ccre_type] += bp_overlap[TEclass][ccre_type]
			if ccre_type not in all_nums:
				all_nums[ccre_type] = 0
				all_bps[ccre_type] = 0
			all_nums[ccre_type] += num_overlap[TEclass][ccre_type]
			all_bps[ccre_type] += bp_overlap[TEclass][ccre_type]
	for ccre_type in TE_nums:
		o.write('TE\t' + ccre_type + '\t' + str(TE_nums[ccre_type]) + '\t' + str(TE_bps[ccre_type]) + '\t' + 
			str(TE_nums[ccre_type]/num_ccres[ccre_type]) + '\t' + str(TE_bps[ccre_type]/bp_ccres[ccre_type]) + '\n')
	for ccre_type in all_nums:
		o.write('All_REs\t' + ccre_type + '\t' + str(all_nums[ccre_type]) + '\t' + str(all_bps[ccre_type]) + '\t' + 
			str(all_nums[ccre_type]/num_ccres[ccre_type]) + '\t' + str(all_bps[ccre_type]/bp_ccres[ccre_type]) + '\n')
	for ccre_type in all_nums:
		o.write('non-TE\t' + ccre_type + '\t' + str(num_ccres[ccre_type] - TE_nums[ccre_type]) + '\t' + str(bp_ccres[ccre_type] - TE_bps[ccre_type]) + '\t' + 
			str(1-(TE_nums[ccre_type]/num_ccres[ccre_type])) + '\t' + str(1-(TE_bps[ccre_type]/bp_ccres[ccre_type])) + '\n')

