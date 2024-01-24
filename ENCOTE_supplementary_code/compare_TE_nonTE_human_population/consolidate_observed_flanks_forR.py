'''
Combines observed and flanks cCRE variant frequency info for R plotting.
Usage: python3 consolidate_observed_flanks_forR.py <observed varFreq file> <left flank varFreq file> <right flank varFreq file> <output file basename>
'''

import sys
import numpy as np

if len(sys.argv) != 5:
	sys.exit(__doc__)

#exclude_top_percentile = 95

zeros = {'observed' : {'TE': {'PLS': [0, 0], 'pELS': [0, 0], 'dELS': [0, 0], 'DNase-H3K4me3': [0, 0], 'CTCF-only': [0, 0]},
		'non-TE': {'PLS': [0, 0], 'pELS': [0, 0], 'dELS': [0, 0], 'DNase-H3K4me3': [0, 0], 'CTCF-only': [0, 0]}},
	'left_flank': {'TE': {'PLS': [0, 0], 'pELS': [0, 0], 'dELS': [0, 0], 'DNase-H3K4me3': [0, 0], 'CTCF-only': [0, 0]},
		'non-TE': {'PLS': [0, 0], 'pELS': [0, 0], 'dELS': [0, 0], 'DNase-H3K4me3': [0, 0], 'CTCF-only': [0, 0]}},
	'right_flank': {'TE': {'PLS': [0, 0], 'pELS': [0, 0], 'dELS': [0, 0], 'DNase-H3K4me3': [0, 0], 'CTCF-only': [0, 0]},
		'non-TE': {'PLS': [0, 0], 'pELS': [0, 0], 'dELS': [0, 0], 'DNase-H3K4me3': [0, 0], 'CTCF-only': [0, 0]}}}

di = {}
observed = {'TE': {'PLS': [], 'pELS': [], 'dELS': [], 'DNase-H3K4me3': [], 'CTCF-only': []},
		'non-TE': {'PLS': [], 'pELS': [], 'dELS': [], 'DNase-H3K4me3': [], 'CTCF-only': []}}
with open(sys.argv[1], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		ccre_id = fields[3]
		ccre = fields[5].replace(',CTCF-bound', '')
		varFreq = float(fields[-1])
		if fields[6] == 'NA':
			if varFreq == 0:
				zeros['observed']['non-TE'][ccre][0] += 1
			else:
				observed['non-TE'][ccre].append(varFreq)
				di[ccre_id] = ['non-TE', ccre, varFreq]
			zeros['observed']['non-TE'][ccre][1] += 1
		else:
			if varFreq == 0:
				zeros['observed']['TE'][ccre][0] += 1
			else:
				observed['TE'][ccre].append(varFreq)
				di[ccre_id] = ['TE', ccre, varFreq]
			zeros['observed']['TE'][ccre][1] += 1
#Write observed cCRE variant frequency for plotting distribution (excluding top percentile to remove long tail)
with open(sys.argv[-1] + '_forDistribution', 'w') as o:
	o.write('Region\tAnnotation\tcCRE_type\tvarFreq\n')
	#thresholds = {}
	#for annotation in observed:
	#	thresholds[annotation] = {}
	#	for ccre in observed[annotation]:
	#		threshold = np.percentile(observed[annotation][ccre], exclude_top_percentile)
	#		thresholds[annotation][ccre] = threshold
	for ccre_id in sorted(di):
		annotation, ccre, varFreq = di[ccre_id]
		o.write('cCRE\t' + annotation + '\t' + ccre + '\t' + str(varFreq) + '\n')
	#	if varFreq < thresholds[annotation][ccre]:
	#		o.write('cCRE\t' + annotation + '\t' + ccre + '\t' + str(varFreq) + '\n')

di = {}
left_flank = {'TE': {'PLS': [], 'pELS': [], 'dELS': [], 'DNase-H3K4me3': [], 'CTCF-only': []},
		'non-TE': {'PLS': [], 'pELS': [], 'dELS': [], 'DNase-H3K4me3': [], 'CTCF-only': []}}
with open(sys.argv[2], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		ccre_id = fields[3]
		ccre = fields[5].replace(',CTCF-bound', '')
		varFreq = float(fields[-1])
		if fields[6] == 'NA':
			if varFreq == 0:
				zeros['left_flank']['non-TE'][ccre][0] += 1
			else:
				left_flank['non-TE'][ccre].append(varFreq)
				di[ccre_id] = ['non-TE', ccre, varFreq]
			zeros['left_flank']['non-TE'][ccre][1] += 1
		else:
			if varFreq == 0:
				zeros['left_flank']['TE'][ccre][0] += 1
			else:
				left_flank['TE'][ccre].append(varFreq)
				di[ccre_id] = ['TE', ccre, varFreq]
			zeros['left_flank']['TE'][ccre][1] += 1
#Write left flank cCRE variant frequency for plotting distribution (excluding top percentile to remove long tail)
with open(sys.argv[-1] + '_forDistribution', 'a') as o:
	#thresholds = {}
	#for annotation in left_flank:
	#	thresholds[annotation] = {}
	#	for ccre in left_flank[annotation]:
	#		threshold = np.percentile(left_flank[annotation][ccre], exclude_top_percentile)
	#		thresholds[annotation][ccre] = threshold
	for ccre_id in sorted(di):
		annotation, ccre, varFreq = di[ccre_id]
		if annotation == 'non-TE':
			continue
		o.write('Left_flank\tflank\t' + ccre + '\t' + str(varFreq) + '\n')
	#	if varFreq < thresholds[annotation][ccre]:
	#		o.write('Left_flank\tflank\t' + ccre + '\t' + str(varFreq) + '\n')

di = {}
right_flank = {'TE': {'PLS': [], 'pELS': [], 'dELS': [], 'DNase-H3K4me3': [], 'CTCF-only': []},
		'non-TE': {'PLS': [], 'pELS': [], 'dELS': [], 'DNase-H3K4me3': [], 'CTCF-only': []}}
with open(sys.argv[3], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		ccre_id = fields[3]
		ccre = fields[5].replace(',CTCF-bound', '')
		varFreq = float(fields[-1])
		if fields[6] == 'NA':
			if varFreq == 0:
				zeros['right_flank']['non-TE'][ccre][0] += 1
			else:
				right_flank['non-TE'][ccre].append(varFreq)
				di[ccre_id] = ['non-TE', ccre, varFreq]
			zeros['right_flank']['non-TE'][ccre][1] += 1
		else:
			if varFreq == 0:
				zeros['right_flank']['TE'][ccre][0] += 1
			else:
				right_flank['TE'][ccre].append(varFreq)
				di[ccre_id] = ['TE', ccre, varFreq]
			zeros['right_flank']['TE'][ccre][1] += 1
#Write right flank cCRE variant frequency for plotting distribution (excluding top percentile to remove long tail)
with open(sys.argv[-1] + '_forDistribution', 'a') as o:
	#thresholds = {}
	#for annotation in right_flank:
	#	thresholds[annotation] = {}
	#	for ccre in right_flank[annotation]:
	#		threshold = np.percentile(right_flank[annotation][ccre], exclude_top_percentile)
	#		thresholds[annotation][ccre] = threshold
	for ccre_id in sorted(di):
		annotation, ccre, varFreq = di[ccre_id]
		if annotation == 'non-TE':
			continue
		o.write('Right_flank\tflank\t' + ccre + '\t' + str(varFreq) + '\n')
	#	if varFreq < thresholds[annotation][ccre]:
	#		o.write('Right_flank\tflank\t' + ccre + '\t' + str(varFreq) + '\n')

#Write observed and shuffled cCRE percent that ovelrap with variants (100% minus percent with varFreq = 0)
with open(sys.argv[-1] + '_forOverlapPercent', 'w') as o:
	o.write('Region\tAnnotation\tcCRE_type\tPercent\n')
	for region in sorted(zeros):
		for annotation in sorted(zeros[region]):
			if "flank" in region:
				if annotation == 'non-TE':
					continue
				else:
					for ccre in sorted(zeros[region][annotation]):
						o.write(region + '\tflank\t' + ccre + '\t' + str((1-zeros[region][annotation][ccre][0]/zeros[region][annotation][ccre][1])*100) + '\n')
			else:
				for ccre in sorted(zeros[region][annotation]):
					o.write(region + '\t' + annotation + '\t' + ccre + '\t' + str((1-zeros[region][annotation][ccre][0]/zeros[region][annotation][ccre][1])*100) + '\n')

