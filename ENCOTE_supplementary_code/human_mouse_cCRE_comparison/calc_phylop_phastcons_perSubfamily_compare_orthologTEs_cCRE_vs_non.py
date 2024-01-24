'''
Compares orthologous TE-cCRE vs. non-orthologous TE-cCRE elements in the same TE subfamily for difference in phyloP and phastCons scores.
Only looks at TEs that are shared between species and are cCREs in one species (anchored from there).
Specifically, uses Mann-Whitney U test to test whether orthologous TE-cCRE and non-orthologous TE-cCRE have the same distributions of phyloP and phastCons scores.
Usage: python3 calc_phylop_phastcons_perSubfamily_compare_orthologTEs_cCRE_vs_non.py <cCRE annotation file> <all TE phylop scores file> <all TE phastcons scores file> <output file>
'''

import sys
import numpy as np
from scipy.stats import mannwhitneyu

if len(sys.argv) != 5:
	sys.exit(__doc__)

#Coverage required for using phyloP or phastCons score
cov_thresh = 0.5 #At least 50% of the TE has scores
#False discovery rate threshold (for determining significance based on q-value)
fdr = 0.05

def p_adjust_bh(p):
	p = np.asfarray(p)
	by_descend = p.argsort()[::-1]
	by_orig = by_descend.argsort()
	steps = float(len(p)) / np.arange(len(p), 0, -1)
	q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
	return q[by_orig]

di = {}
di_info = {}
with open(sys.argv[1], 'r') as f:
	header = f.readline().strip('\n').split('\t')
	TE_coord_index = header.index('Human_TE')
	TE_subfamily_index = header.index('Human_repName')
	TE_species_index = header.index('Species_classification')
	TE_clade_index = header.index('Clade')
	synteny_annotation_index = header.index('Synteny_annotation')
	ccre_comparison_index = header.index('cCRE_comparison')
	for line in f:
		fields = line.rstrip('\n').split('\t')
		coord = fields[TE_coord_index]
		subfamily = fields[TE_subfamily_index]
		species_classification = fields[TE_species_index]
		clade = fields[TE_clade_index]
		if species_classification != 'Shared_subfamily': #Skip lineage specific subfamilies
			continue
		if subfamily not in di:
			di[subfamily] = {}
			di_info[subfamily] = {}
		if coord not in di[subfamily]:
			di[subfamily][coord] = {}
		di[subfamily][coord]['type'] = 'cCRE'
		di_info[subfamily]['species'] = species_classification
		di_info[subfamily]['clade'] = clade
		di[subfamily][coord]['synteny'] = fields[synteny_annotation_index]
		if fields[ccre_comparison_index] == 'Shared_cCRE' or fields[ccre_comparison_index] == 'Different_cCRE':
			di[subfamily][coord]['ccre'] = 'both_species'
		else:
			di[subfamily][coord]['ccre'] = 'one_species'
#Add phyloP and phastCons scores
with open(sys.argv[2], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		coord = fields[0] + ':' + fields[1] + '-' + fields[2] + '(' + fields[5] + ')'
		subfamily = fields[3]
		if subfamily not in di: #No cCRE elements in the TE subfamily, skip
			continue
		if coord not in di[subfamily]:
			continue #Skip if not cCRE
#			di[subfamily][coord] = {}
#			di[subfamily][coord]['type'] = 'non-cCRE'
		if fields[6] == 'NA':
			di[subfamily][coord]['phylop_score'] = 'NA'
		else:
			di[subfamily][coord]['phylop_score'] = float(fields[6])
		di[subfamily][coord]['phylop_bases'] = int(fields[7])
		di[subfamily][coord]['phylop_cov'] = float(fields[8])
with open(sys.argv[3], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		coord = fields[0] + ':' + fields[1] + '-' + fields[2] + '(' + fields[5] + ')'
		subfamily = fields[3]
		if subfamily not in di: #No cCRE elements in the TE subfamily, skip
			continue
		if coord not in di[subfamily]:
			continue #Skip if not cCRE
#			di[subfamily][coord] = {}
#			di[subfamily][coord]['type'] = 'non-cCRE'
		if fields[6] == 'NA':
			di[subfamily][coord]['phastcons_score'] = 'NA'
		else:
			di[subfamily][coord]['phastcons_score'] = float(fields[6])
		di[subfamily][coord]['phastcons_bases'] = int(fields[7])
		di[subfamily][coord]['phastcons_cov'] = float(fields[8])

#Perform Mann-Whitney U test
#Print all classifications to standard output (can be used to do Mann-Whitney U test and/or plot across all TE subfamilies)
phylop_pvals = []
phastcons_pvals = []
list_subfamilies = sorted(di)
list_tests = {'phylop': [], 'phastcons': []}
di_quartiles = {}
di_numbers = {}
print('Subfamily\tScore_type\tGroup\tScore\tClade')
for subfamily in list_subfamilies:
	#Group12 is orthologous TE + cCRE in human and mouse
	#Group3 is non-orthologous TE (cCRE in one species)
	group12_phylop = []
	group12_phastcons = []
	group3_phylop = []
	group3_phastcons = []
	for coord in di[subfamily]:
		if di[subfamily][coord]['synteny'] == 'TEortholog':
			if di[subfamily][coord]['ccre'] == 'both_species':
				if 'phylop_score' in di[subfamily][coord]:
					if di[subfamily][coord]['phylop_cov'] >= cov_thresh:
						group12_phylop.append(di[subfamily][coord]['phylop_score'])
						print(subfamily + '\tphyloP\tgroup12\t' + str(di[subfamily][coord]['phylop_score']) + '\t' + di_info[subfamily]['clade'])
				if 'phastcons_score' in di[subfamily][coord]:
					if di[subfamily][coord]['phastcons_cov'] >= cov_thresh:
						group12_phastcons.append(di[subfamily][coord]['phastcons_score'])
						print(subfamily + '\tphastCons\tgroup12\t' + str(di[subfamily][coord]['phastcons_score']) + '\t' + di_info[subfamily]['clade'])
			else:
				if 'phylop_score' in di[subfamily][coord]:
					if di[subfamily][coord]['phylop_cov'] >= cov_thresh:
						group12_phylop.append(di[subfamily][coord]['phylop_score'])
						print(subfamily + '\tphyloP\tgroup12\t' + str(di[subfamily][coord]['phylop_score']) + '\t' + di_info[subfamily]['clade'])
				if 'phastcons_score' in di[subfamily][coord]:
					if di[subfamily][coord]['phastcons_cov'] >= cov_thresh:
						group12_phastcons.append(di[subfamily][coord]['phastcons_score'])
						print(subfamily + '\tphastCons\tgroup12\t' + str(di[subfamily][coord]['phastcons_score']) + '\t' + di_info[subfamily]['clade'])
		elif di[subfamily][coord]['synteny'] == 'Human_TE':
			if di[subfamily][coord]['ccre'] == 'one_species':
				if 'phylop_score' in di[subfamily][coord]:
					if di[subfamily][coord]['phylop_cov'] >= cov_thresh:
						group3_phylop.append(di[subfamily][coord]['phylop_score'])
						print(subfamily + '\tphyloP\tgroup3\t' + str(di[subfamily][coord]['phylop_score']) + '\t' + di_info[subfamily]['clade'])
				if 'phastcons_score' in di[subfamily][coord]:
					if di[subfamily][coord]['phastcons_cov'] >= cov_thresh:
						group3_phastcons.append(di[subfamily][coord]['phastcons_score'])
						print(subfamily + '\tphastCons\tgroup3\t' + str(di[subfamily][coord]['phastcons_score']) + '\t' + di_info[subfamily]['clade'])
	#See if there are enough elements to do the test for phyloP
	if len(group12_phylop) < 5 or len(group3_phylop) < 5:
		continue #Fewer than 5 elements in any group, skip
	di_quartiles[subfamily] = {}
	di_numbers[subfamily] = {}
	#phyloP Mann-Whitney U test
	ustat, pval = mannwhitneyu(group12_phylop, group3_phylop, alternative='two-sided')
	phylop_pvals.append(pval)
	list_tests['phylop'].append(subfamily + '/group12/group3')
	di_quartiles[subfamily]['phylop'] = {}
	di_quartiles[subfamily]['phylop']['group12'] = np.quantile(group12_phylop, [0.25, 0.5, 0.75])
	di_quartiles[subfamily]['phylop']['group3'] = np.quantile(group3_phylop, [0.25, 0.5, 0.75])
	di_numbers[subfamily]['phylop'] = {}
	di_numbers[subfamily]['phylop']['group12'] = len(group12_phylop)
	di_numbers[subfamily]['phylop']['group3'] = len(group3_phylop)
	#See if there are enough elements to do the test for phastCons
	if len(group12_phastcons) < 5 or len(group3_phastcons) < 5:
		continue #Fewer than 5 elements in any group, skip
	#phastCons Mann-Whitney U test
	ustat, pval = mannwhitneyu(group12_phastcons, group3_phastcons, alternative='two-sided')
	phastcons_pvals.append(pval)
	list_tests['phastcons'].append(subfamily + '/group12/group3')
	di_quartiles[subfamily]['phastcons'] = {}
	di_quartiles[subfamily]['phastcons']['group12'] = np.quantile(group12_phastcons, [0.25, 0.5, 0.75])
	di_quartiles[subfamily]['phastcons']['group3'] = np.quantile(group3_phastcons, [0.25, 0.5, 0.75])
	di_numbers[subfamily]['phastcons'] = {}
	di_numbers[subfamily]['phastcons']['group12'] = len(group12_phastcons)
	di_numbers[subfamily]['phastcons']['group3'] = len(group3_phastcons)
#Benjamini-Hochberg correction
phylop_qvals = p_adjust_bh(phylop_pvals)
phastcons_qvals = p_adjust_bh(phastcons_pvals)

quantiles = [str(25), str(50), str(75)]
groups = ['group12', 'group3']
summary_di = {}
#Write to output
with open(sys.argv[-1], 'w') as o:
	o.write('Subfamily\tScore_type\tGroup12_number\tGroup3_number\tCall\tq-value\tp-value\tMean_difference\tClade\n')
	#phyloP test result
	for i in range(len(list_tests['phylop'])):
		subfamily, groupA, groupB = list_tests['phylop'][i].split('/')
		if subfamily not in summary_di:
			summary_di[subfamily] = {}
			summary_di[subfamily]['phylop'] = {}
		comparison = groupA + '/' + groupB
		summary_di[subfamily]['phylop'][comparison] = {}
		phylop_qval = phylop_qvals[i]
		phylop_pval = phylop_pvals[i]
		phylop_diff = []
		for j in range(3):
			phylop_diff.append(di_quartiles[subfamily]['phylop'][groupA][j] - di_quartiles[subfamily]['phylop'][groupB][j])
		mean_phylop_diff = sum(phylop_diff)/len(phylop_diff)
		summary_di[subfamily]['phylop'][comparison]['phylop_diff'] = mean_phylop_diff
		o.write(subfamily + '\tphyloP\t' + str(di_numbers[subfamily]['phylop'][groupA]) + '\t' + str(di_numbers[subfamily]['phylop'][groupB]) + '\t')
		if phylop_qval < fdr:
			if mean_phylop_diff == 'NA':
				o.write('NA')
				summary_di[subfamily]['phylop'][comparison]['call_phylop'] = 'NA'
			elif mean_phylop_diff > 0:
				o.write('Greater')
				summary_di[subfamily]['phylop'][comparison]['call_phylop'] = 'Greater'
			else:
				o.write('Less')
				summary_di[subfamily]['phylop'][comparison]['call_phylop'] = 'Less'
		else:
			o.write('Not_sig')
			summary_di[subfamily]['phylop'][comparison]['call_phylop'] = 'Not_sig'
		o.write('\t' + str(phylop_qval) + '\t' + str(phylop_pval) + '\t' + str(mean_phylop_diff) + '\t' + di_info[subfamily]['clade'] + '\n')
	#phastCons test results
	for i in range(len(list_tests['phastcons'])):
		subfamily, groupA, groupB = list_tests['phastcons'][i].split('/')
		if subfamily not in summary_di:
			summary_di[subfamily] = {}
		if 'phastcons' not in summary_di[subfamily]:
			summary_di[subfamily]['phastcons'] = {}
		comparison = groupA + '/' + groupB
		summary_di[subfamily]['phastcons'][comparison] = {}
		phastcons_qval = phastcons_qvals[i]
		phastcons_pval = phastcons_pvals[i]
		phastcons_diff = []
		for j in range(3):
			phastcons_diff.append(di_quartiles[subfamily]['phastcons'][groupA][j] - di_quartiles[subfamily]['phastcons'][groupB][j])
		mean_phastcons_diff = sum(phastcons_diff)/len(phastcons_diff)
		summary_di[subfamily]['phastcons'][comparison]['phastcons_diff'] = mean_phastcons_diff
		o.write(subfamily + '\tphastCons\t' + str(di_numbers[subfamily]['phastcons'][groupA]) + '\t' + str(di_numbers[subfamily]['phastcons'][groupB]) + '\t')
		if phastcons_qval < fdr:
			if mean_phastcons_diff == 'NA':
				o.write('NA')
				summary_di[subfamily]['phastcons'][comparison]['call_phastcons'] = 'NA'
			elif mean_phastcons_diff > 0:
				o.write('Greater')
				summary_di[subfamily]['phastcons'][comparison]['call_phastcons'] = 'Greater'
			else:
				o.write('Less')
				summary_di[subfamily]['phastcons'][comparison]['call_phastcons'] = 'Less'
		else:
			o.write('Not_sig')
			summary_di[subfamily]['phastcons'][comparison]['call_phastcons'] = 'Not_sig'
		o.write('\t' + str(phastcons_qval) + '\t' + str(phastcons_pval) + '\t' + str(mean_phastcons_diff) + '\t' + di_info[subfamily]['clade'] + '\n')

with open(sys.argv[-1] + '_quartiles', 'w') as o:
	o.write('Subfamily\tScore_type\tGroup\tQuantile\tScore\tClade\n')
	for subfamily in list_subfamilies:
		if subfamily in di_quartiles:
			if 'phylop' in di_quartiles[subfamily]:
				for group in groups:
					for j in range(len(quantiles)):
						o.write(subfamily + '\tphyloP\t' + group + '\t' + quantiles[j] + '\t' + str(di_quartiles[subfamily]['phylop'][group][j]) + 
							'\t' + di_info[subfamily]['clade'] + '\n')
			if 'phastcons' in di_quartiles[subfamily]:
				for group in groups:
					for j in range(len(quantiles)):
						o.write(subfamily + '\tphastCons\t' + group + '\t' + quantiles[j] + '\t' + str(di_quartiles[subfamily]['phastcons'][group][j]) + 
							'\t' + di_info[subfamily]['clade'] + '\n')

#Print summary of each subfamily to standard output
#print('Subfamily\tScore_type\tCall\tMean_quantile_diff\tSpecies_classification\tClade')
#for subfamily in list_subfamilies:
#	print(subfamily + '\tphyloP\t' + summary_di[subfamily]['call_phylop'] + '\t' + str(summary_di[subfamily]['phylop']) + '\t' + 
#		di_info[subfamily]['species'] + '\t' + di_info[subfamily]['clade'])
#	print(subfamily + '\tphastCons\t' + summary_di[subfamily]['call_phastcons'] + '\t' + str(summary_di[subfamily]['phastcons']) + '\t' + 
#		di_info[subfamily]['species'] + '\t' + di_info[subfamily]['clade'])

