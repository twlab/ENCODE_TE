'''
Summarizes results from GWAS-cCRE overlap.
<GWAS traits file> used to group disease/traits based on parent term (e.g. cancer, immune, etc)
This version is for shuffled cCRE coordinates (used as random background). Only difference is that there is standard output.
Usage: python3 summarize_gwas_ccre_overlap_shuffled.py <input overlap file> <GWAS traits file> <output file>
'''

import sys

if len(sys.argv) != 4:
	sys.exit(__doc__)

counts_di = {'all': {'ccre': 0, 'nonTE': 0, 'TE': 0, 'human_nonTE': 0, 'human_TE': 0}}
trait_di = {}
efo_di = {'uri': {}, 'name': {}}
with open(sys.argv[2], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		disease_trait = fields[0]
		efo_term = fields[1]
		efo_uri = fields[2]
		parent_term = fields[3]
		parent_uri = fields[4]
		if parent_term not in counts_di:
			counts_di[parent_term] = {'ccre': 0, 'nonTE': 0, 'TE': 0, 'human_nonTE': 0, 'human_TE': 0}
		if disease_trait not in trait_di:
			trait_di[disease_trait] = [parent_term]
		else:
			trait_di[disease_trait].append(parent_term)
		if efo_term not in efo_di['name']:
			efo_di['name'][efo_term] = [parent_term]
		else:
			efo_di['name'][efo_term].append(parent_term)
		if efo_uri not in efo_di['uri']:
			efo_di['uri'][efo_uri] = [parent_term]
		else:
			efo_di['uri'][efo_uri].append(parent_term)

already_counted = {}
for parent_term in counts_di:
	already_counted[parent_term] = {}
with open(sys.argv[1], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		coords = fields[0] + ':' + fields[1] + '-' + fields[2]
		disease_trait = fields[3]
		efo_terms = fields[5].split(',')
		efo_uris = fields[6].split(',')
		list_parent_terms = []
		if disease_trait in trait_di:
			list_parent_terms = trait_di[disease_trait]
		else:
			for i in range(len(efo_uris)):
				efo_uri = efo_uris[i].strip()
				efo_term = efo_terms[i].strip()
				if efo_uri in efo_di['uri']:
					parent_term = efo_di['uri'][efo_uri][i]
					list_parent_terms.append(parent_term)
				elif efo_term in efo_di['name']:
					parent_term = efo_di['name'][efo_term][i]
					list_parent_terms.append(parent_term)
				else:
					continue
		terms_to_add = []
		not_added = True
		if coords not in already_counted:
			already_counted[coords] = {'all': 'Y'}
			for parent_term in set(list_parent_terms):
				already_counted[coords][parent_term] = 'Y'
				terms_to_add.append(parent_term)
		else:
			for parent_term in set(list_parent_terms):
				if parent_term in already_counted[coords]: #Variant already counted once for current parent term, skip
					continue
				else:
					terms_to_add.append(parent_term)
					already_counted[coords][parent_term] = 'Y'
			not_added = False
		if terms_to_add == []:
			continue
		ccre_chr = fields[7].strip()
		TE_chr = fields[13].strip()
		synteny_info = fields[19]
		if not_added:
			if ccre_chr != '.':
				counts_di['all']['ccre'] += 1
				if TE_chr == 'NA':
					counts_di['all']['nonTE'] += 1
				else:
					counts_di['all']['TE'] += 1
				if synteny_info == 'Human_specific':
					counts_di['all']['human_nonTE'] += 1
				elif synteny_info == 'Human_TE':
					counts_di['all']['human_TE'] += 1
		for parent_term in terms_to_add:
			if parent_term not in counts_di:
				sys.stderr.write(parent_term + ' not in counts_di for some reason\n')
				sys.exit()
			if ccre_chr != '.':
				counts_di[parent_term]['ccre'] += 1
			else:
				continue
			if TE_chr == 'NA':
				counts_di[parent_term]['nonTE'] += 1
			else:
				counts_di[parent_term]['TE'] += 1
			if synteny_info == 'Human_specific':
				counts_di[parent_term]['human_nonTE'] += 1
			elif synteny_info == 'Human_TE':
				counts_di[parent_term]['human_TE'] += 1

with open(sys.argv[-1], 'w') as o:
	o.write('Parent_term\tcCRE_counts\tnonTE_counts\tTE_counts\tHuman_nonTE_counts\tHuman_TE_counts\n')
	o.write('All\t' + str(counts_di['all']['ccre']) + '\t' + str(counts_di['all']['nonTE']) + '\t' + str(counts_di['all']['TE']) + '\t' + 
		str(counts_di['all']['human_nonTE']) + '\t' + str(counts_di['all']['human_TE']) + '\n')
	for parent_term in sorted(counts_di):
		if parent_term != 'all':
			o.write(parent_term + '\t' + str(counts_di[parent_term]['ccre']) + '\t' + str(counts_di[parent_term]['nonTE']) + '\t' + str(counts_di[parent_term]['TE']) + '\t' +
				str(counts_di[parent_term]['human_nonTE']) + '\t' + str(counts_di[parent_term]['human_TE']) + '\n')

