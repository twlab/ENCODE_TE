'''
Call likely turnover events using phastCons scores of top nearby TF binding loci.
Uses mean/median phastCons for syntenic TF binding regions to infer which TF binding sites are ancestral.
	If either the reference or syntenic site are TE-associated, use lineage-specific TE annotation to improve ancestry inference.
	This version (v2) uses orthology/synteny of TF binding sites to infer ancestry (doesn't actually use reference/syntenic lineage TE info).
Usage: python3 call_likely_turnover_v2.py <top binding site loci with phastCons> <syntenic TF binding phastCons summary> <reference lineage TEs> <syntenic lineage TEs> <output file basename>
'''

import sys

if len(sys.argv) != 6:
	sys.exit(__doc__)

ref_TEs = {}
with open(sys.argv[3], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		clade = fields[1]
		ref_TEs[subfamily] = clade
syn_TEs = {}
with open(sys.argv[4], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		clade = fields[1]
		syn_TEs[subfamily] = clade

syntenic_thresholds = {}
with open(sys.argv[2], 'r') as f:
	header = f.readline().rstrip('\n').split('\t')
	mean_index = header.index('Mean')
	median_index = header.index('Median')
	for line in f:
		fields = line.rstrip('\n').split('\t')
		TF = fields[0]
		genome = fields[1]
		mean = float(fields[mean_index])
		median = float(fields[median_index])
		if TF not in syntenic_thresholds:
			syntenic_thresholds[TF] = {}
		syntenic_thresholds[TF][genome] = {'mean': mean, 'median': median}

list_ref_peaks = []
turnover_di = {'mean': {}, 'median': {}, 'summary': {}}
with open(sys.argv[1], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		if fields[-2] == 'NA' or fields[-1] == 'NA':
			continue
		ref_peak = fields[4]
		list_ref_peaks.append(ref_peak)
		ref_peak_orthology = ref_peak.split('|')[3]
		TF = fields[-3]
		if TF not in turnover_di['summary']:
			turnover_di['summary'][TF] = {'mean': {'neither': 0, 'ref': 0, 'syn': 0, 'both': 0}, 'median': {'neither': 0, 'ref': 0, 'syn': 0, 'both': 0}}
		ref_phastcons = float(fields[-2])
		syn_peaks = fields[5].split(',')
		syn_phastcons = max([float(x) for x in fields[-1].split(',')])
		syn_index = fields[-1].split(',').index(str(syn_phastcons))
		syn_peak = syn_peaks[syn_index]
		syn_peak_orthology = syn_peak.split('|')[3]
		for sum_type in syntenic_thresholds[TF]['Reference']:
			if ref_peak_orthology == 'Orth' and syn_peak_orthology == 'Orth': #Both reference and syntenic sites have orthology/synteny with the other species
				if ref_phastcons < syntenic_thresholds[TF]['Reference'][sum_type] and syn_phastcons < syntenic_thresholds[TF]['Syntenic'][sum_type]:
					turnover_di[sum_type][ref_peak] = line.rstrip('\n') + '\tNeither\n'
					turnover_di['summary'][TF][sum_type]['neither'] += 1
				elif ref_phastcons >= syntenic_thresholds[TF]['Reference'][sum_type] and syn_phastcons >= syntenic_thresholds[TF]['Syntenic'][sum_type]:
					turnover_di[sum_type][ref_peak] = line.rstrip('\n') + '\tBoth\n'
					turnover_di['summary'][TF][sum_type]['both'] += 1
				elif ref_phastcons >= syntenic_thresholds[TF]['Reference'][sum_type]:
					turnover_di[sum_type][ref_peak] = line.rstrip('\n') + '\tReference\n'
					turnover_di['summary'][TF][sum_type]['ref'] += 1
				elif syn_phastcons >= syntenic_thresholds[TF]['Syntenic'][sum_type]:
					turnover_di[sum_type][ref_peak] = line.rstrip('\n') + '\tSyntenic\n'
					turnover_di['summary'][TF][sum_type]['syn'] += 1
			elif ref_peak_orthology == 'Orth': #Only reference site has orthology/synteny => only reference site can be ancestral
				if ref_phastcons >= syntenic_thresholds[TF]['Reference'][sum_type]:
					turnover_di[sum_type][ref_peak] = line.rstrip('\n') + '\tReference\n'
					turnover_di['summary'][TF][sum_type]['ref'] += 1
				else:
					turnover_di[sum_type][ref_peak] = line.rstrip('\n') + '\tNeither\n'
					turnover_di['summary'][TF][sum_type]['neither'] += 1
			elif syn_peak_orthology == 'Orth': #Only syntenic site has orthology/synteny => only syntenic site can be ancestral
				if syn_phastcons >= syntenic_thresholds[TF]['Syntenic'][sum_type]:
					turnover_di[sum_type][ref_peak] = line.rstrip('\n') + '\tSyntenic\n'
					turnover_di['summary'][TF][sum_type]['syn'] += 1
				else:
					turnover_di[sum_type][ref_peak] = line.rstrip('\n') + '\tNeither\n'
					turnover_di['summary'][TF][sum_type]['neither'] += 1
			else: #Neither reference nor syntenic sites have orthology/synteny with the other species => both are lineage specific
				turnover_di[sum_type][ref_peak] = line.rstrip('\n') + '\tNeither\n'
				turnover_di['summary'][TF][sum_type]['neither'] += 1

with open(sys.argv[-1] + '_mean', 'w') as o:
	o.write(header.rstrip('\n') + '\tAncestral_call\n')
	for ref_peak in list_ref_peaks:
		o.write(turnover_di['mean'][ref_peak])
with open(sys.argv[-1] + '_median', 'w') as o:
	o.write(header.rstrip('\n') + '\tAncestral_call\n')
	for ref_peak in list_ref_peaks:
		o.write(turnover_di['median'][ref_peak])
with open(sys.argv[-1] + '_summary', 'w') as o:
	o.write('TF\tThreshold_type\tRef_perc\tRef_num\tSyn_perc\tSyn_num\tBoth_perc\tBoth_num\tNeither_perc\tNeither_num\n')
	for TF in sorted(turnover_di['summary']):
		for sum_type in sorted(turnover_di['summary'][TF]):
			total = (turnover_di['summary'][TF][sum_type]['neither'] + turnover_di['summary'][TF][sum_type]['ref'] + 
				turnover_di['summary'][TF][sum_type]['syn'] + turnover_di['summary'][TF][sum_type]['both'])
			line = (TF + '\t' + sum_type + '\t' + str(round(turnover_di['summary'][TF][sum_type]['ref']/total*100, 1)) + '\t' + str(turnover_di['summary'][TF][sum_type]['ref']) + '\t' +
				str(round(turnover_di['summary'][TF][sum_type]['syn']/total*100, 1)) + '\t' + str(turnover_di['summary'][TF][sum_type]['syn']) + '\t' +
				str(round(turnover_di['summary'][TF][sum_type]['both']/total*100, 1)) + '\t' + str(turnover_di['summary'][TF][sum_type]['both']) + '\t' +
				str(round(turnover_di['summary'][TF][sum_type]['neither']/total*100, 1)) + '\t' + str(turnover_di['summary'][TF][sum_type]['neither']) + '\n')
			o.write(line)

