'''
Parses the output of the liftOver match output to make sense.
This version is for hg38 cCREs (with TE overlap info) lifted to mouse and then overlapped with mouse TEs.
Usage: python3 parse_liftOver_match_TE.py <input file> <output file>
'''

import sys

if len(sys.argv) != 3:
	sys.exit(__doc__)

with open(sys.argv[1], 'r') as f:
	with open(sys.argv[-1], 'w') as o:
		#Header
		o.write('hg38_cCRE\thg38_id1\thg38_id2\thg38_cCRE_type\thg38_TE\thg38_repName\thg38_repClass/Family\thg38_alignment\thg38_TE_ovelap\tmm10_coord\tmm10_TE\tmm10_repName\tmm10_repClass/Family\tmm10_TE_overlap\n')
		for line in f:
			fields = line.rstrip('\n').split('\t')
			new_line = ''
			#Original cCRE info
			orig_ccre = fields[0] + ':' + fields[1] + '-' + fields[2]
			new_line += orig_ccre + '\t' + '\t'.join(fields[3:6])
			#TE overlap info (if any)
			if fields[6] == '.': #No TE overlap
				orig_te = 'NA'
				new_line += '\t' + orig_te + '\tNA\tNA\tNA\t0'
			else:
				orig_te = fields[6] + ':' + fields[7] + '-' + fields[8] + '(' + fields[11] + ')'
				new_line += '\t' + orig_te + '\t' + '\t'.join(fields[9:11]) + '\t' + '\t'.join(fields[12:14])
			#liftOver coordinates (if any)
			if fields[14] == '.': #No liftOver
					liftover_coords = 'NA'
			else:
				liftover_coords = fields[14] + ':' + fields[15] + '-' + fields[16]
			new_line += '\t' + liftover_coords
			#liftOver genome TE info (if any)
			if fields[17] == '.': #No TE overlap
				liftover_te = 'NA'
				new_line += '\t' + liftover_te + '\tNA\tNA\t0'
			else:
				liftover_te = fields[17] + ':' + fields[18] + '-' + fields[19] + '(' + fields[22] + ')'
				new_line += '\t' + liftover_te + '\t' + '\t'.join(fields[20:22]) + '\t' + fields[-1]
			o.write(new_line + '\n')

