'''
Takes matched human cCRE to mouse cCRE/TE files to quantify syntenic elements.
Adds age of TE subfamily.
Writes mouse-human cCRE-TE info to output file.
Prints summary/statistics to standard output.
This version (v2) considers syntenic TEs in the same TE class/family as orthologous.
This version (mouse_human) uses mouse cCREs as reference for comparison to human.
Usage: python3 quantify_mouse_human_ccre_TEs.py <matched cCRE file> <matched TE file> <kimura divergences human> <kimura divergences mouse> <rmsk bed file with TE family>  <output file>
'''

import sys

if len(sys.argv) != 7:
	sys.exit(__doc__)

human_ages = {}
with open(sys.argv[3], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		age = fields[-1]
		human_ages[subfamily] = age

mouse_ages = {}
with open(sys.argv[4], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		age = fields[-1]
		mouse_ages[subfamily] = age

list_mm10_ccres = []
di = {}
with open(sys.argv[1], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		mm10_ccre = fields[0]
		list_mm10_ccres.append(mm10_ccre)
		di[mm10_ccre] = fields[:9]
		if fields[4] != "NA":
			subfamily = fields[5]
			if subfamily in mouse_ages:
				di[mm10_ccre].append(mouse_ages[subfamily])
			else:
				di[mm10_ccre].append('NA')
		else:
			di[mm10_ccre].append('NA')
		di[mm10_ccre] += fields[9:]

with open(sys.argv[2], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		mm10_ccre = fields[0]
		di[mm10_ccre] += fields[10:]
		if fields[10] != "NA":
			subfamily = fields[11]
			if subfamily in human_ages:
				di[mm10_ccre].append(human_ages[subfamily])
			else:
				di[mm10_ccre].append('NA')
		else:
			di[mm10_ccre].append('NA')

#Criteria for determining if subfamily consensus sequences are sufficiently related based on class
subfamily_toClass = {}
with open(sys.argv[5], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[3]
		TE_class = fields[4]
		if subfamily not in subfamily_toClass:
			subfamily_toClass[subfamily] = TE_class

mouse_cCREs = {}
mouse_TEs = {}
syntenic_regions = []
human_cCREs = {} #must be subset of syntenic regions
human_TEs = {} #must be subset of syntenic regions
human_ccre_types = {}
mouse_ccre_types = {}
for mm10_ccre in list_mm10_ccres:
	fields = di[mm10_ccre]
	#Mouse cCRE info
	mouse_ccre_type = []
	if 'dELS' in fields[3]:
		mouse_ccre_type.append('dELS')
	if 'pELS' in fields[3]:
		mouse_ccre_type.append('pELS')
	if 'PLS' in fields[3]:
		mouse_ccre_type.append('PLS')
	if 'DNase-H3K4me3' in fields[3]:
		mouse_ccre_type.append('DNase-H3K4me3')
	if 'CTCF-only' in fields[3]:
		mouse_ccre_type.append('CTCF-only')
	if 'CTCF' in fields[3]:
		mouse_ccre_type.append('CTCF-bound')
	if 'DNase-only' in fields[3]:
		mouse_ccre_type.append('DNase-only')
	for ccre_type in mouse_ccre_type:
		if ccre_type not in mouse_ccre_types:
			mouse_ccre_types[ccre_type] = []
		mouse_ccre_types[ccre_type].append(mm10_ccre)
	mouse_cCREs[mm10_ccre] = mouse_ccre_type
	#Mouse TE, if any
	if fields[4] != 'NA':
		subfamily = fields[5]
		TE_class = fields[6].split('/')[0]
		try:
			family = fields[6].split('/')[1]
		except IndexError:
			family = 'NA'
		age = fields[9]
		mouse_TEs[mm10_ccre] = [subfamily, family, TE_class, age]
	#Syntenic human region, if any
	if fields[10] != 'NA':
		syntenic_regions.append(mm10_ccre)
	#Syntenic human cCRE, if any
	if fields[14] != 'NA':
		human_ccre_type = []
		if 'dELS' in fields[14]:
			human_ccre_type.append('dELS')
		if 'pELS' in fields[14]:
			human_ccre_type.append('pELS')
		if 'PLS' in fields[14]:
			human_ccre_type.append('PLS')
		if 'DNase-H3K4me3' in fields[14]:
			human_ccre_type.append('DNase-H3K4me3')
		if 'CTCF-only' in fields[14]:
			human_ccre_type.append('CTCF-only')
		if 'CTCF' in fields[14]:
			human_ccre_type.append('CTCF-bound')
		if 'DNase-only' in fields[14]:
			human_ccre_type.append('DNase-only')
		for ccre_type in human_ccre_type:
			if ccre_type not in human_ccre_types:
				human_ccre_types[ccre_type] = []
			human_ccre_types[ccre_type].append(mm10_ccre)
		human_cCREs[mm10_ccre] = human_ccre_type
	#Syntenic human TE, if any
	if fields[16] != 'NA':
		subfamily = fields[17]
		TE_class = fields[18].split('/')[0]
		try:
			family = fields[18].split('/')[1]
		except IndexError:
			family = 'NA'
		age = fields[-1]
		human_TEs[mm10_ccre] = [subfamily, family, TE_class, age]

orthologous_TE_numbers = [0, 0, 0]
orthologous_TEs = {'human': [], 'mouse': [], 'human_class': [], 'mouse_class': []}
syntenic_nonOrtholog_TE = {'human': [], 'mouse': []}
orthologous_TE_ccres = {'total': {'PLS': 0, 'pELS': 0, 'dELS': 0, 'DNase-H3K4me3': 0, 'CTCF-only': 0, 'DNase-only': 0, 'all': 0},
			'same': {'PLS': 0, 'pELS': 0, 'dELS': 0, 'DNase-H3K4me3': 0, 'CTCF-only': 0, 'DNase-only': 0, 'all': 0}, 
			'different': {'PLS': 0, 'pELS': 0, 'dELS': 0, 'DNase-H3K4me3': 0, 'CTCF-only': 0, 'DNase-only': 0, 'all': 0}, 
			'mouse_only': {'PLS': 0, 'pELS': 0, 'dELS': 0, 'DNase-H3K4me3': 0, 'CTCF-only': 0, 'DNase-only': 0, 'all': 0}}
syntenic_nonTE_numbers = {'total': {'PLS': 0, 'pELS': 0, 'dELS': 0, 'DNase-H3K4me3': 0, 'CTCF-only': 0, 'DNase-only': 0, 'all': 0},
				'same': {'PLS': 0, 'pELS': 0, 'dELS': 0, 'DNase-H3K4me3': 0, 'CTCF-only': 0, 'DNase-only': 0, 'all': 0}, 
				'different': {'PLS': 0, 'pELS': 0, 'dELS': 0, 'DNase-H3K4me3': 0, 'CTCF-only': 0, 'DNase-only': 0, 'all': 0}, 
				'mouse_only': {'PLS': 0, 'pELS': 0, 'dELS': 0, 'DNase-H3K4me3': 0, 'CTCF-only': 0, 'DNase-only': 0, 'all': 0}}
syntenic_regions_set = set(syntenic_regions)
#For species-specific cCREs
list_nonHuman_TEs = []
novel_cCREs = {'TE': {'PLS': 0, 'pELS': 0, 'dELS': 0, 'DNase-H3K4me3': 0, 'CTCF-only': 0, 'DNase-only': 0, 'all': 0},
		'non-TE': {'PLS': 0, 'pELS': 0, 'dELS': 0, 'DNase-H3K4me3': 0, 'CTCF-only': 0, 'DNase-only': 0, 'all': 0}}
shared_cCREs = {'TE': {'PLS': 0, 'pELS': 0, 'dELS': 0, 'DNase-H3K4me3': 0, 'CTCF-only': 0, 'DNase-only': 0, 'all': 0},
		'non-TE': {'PLS': 0, 'pELS': 0, 'dELS': 0, 'DNase-H3K4me3': 0, 'CTCF-only': 0, 'DNase-only': 0, 'all': 0}}
#Write info to output file
with open(sys.argv[-1], 'w') as o:
	o.write('Mouse_cCRE\tMouse_cCRE_id1\tMouse_cCRE_id2\tMouse_cCRE_type\tMouse_TE\tMouse_repName\tMouse_repClass/Family\talignment\tMouse_TE_overlap\tMouse_TE_age\t' +
		'Human_syntenic\tHuman_cCRE\tHuman_cCRE_id1\tHuman_cCRE_id2\tHuman_cCRE_type\tHuman_cCRE_overlap\tHuman_TE\tHuman_repName\tHuman_repClass/Family\tHuman_TE_overlap\tHuman_TE_age\t' +
		'Synteny_annotation\tcCRE_comparison\n')
	for mm10_ccre in list_mm10_ccres:
		o.write('\t'.join(di[mm10_ccre]))
		#Quantify orthologous TE (shared TE) contribution to conserved/different/species-specific (novel) cCREs
		#Get numbers for orthologous TEs between human and mouse in the meantime
		if mm10_ccre in human_TEs and mm10_ccre in mouse_TEs:
			orthologous_TE_numbers[2] += 1
			mouse_TE = di[mm10_ccre][4]
			human_TE = di[mm10_ccre][16]
			mouse_subfamily = human_TEs[mm10_ccre][0]
			human_subfamily = mouse_TEs[mm10_ccre][0]
			mouse_age = human_TEs[mm10_ccre][3]
			human_age = mouse_TEs[mm10_ccre][3]
			is_ortholog = False
			if mouse_subfamily == human_subfamily:
				orthologous_TE_numbers[0] += 1
				orthologous_TEs['human'].append([human_TE, human_subfamily, str(human_age)])
				orthologous_TEs['mouse'].append([mouse_TE, mouse_subfamily, str(mouse_age)])
				is_ortholog = True
			#Consider subfamilies that share the same TE class as the same (obviously only if the subfamily is shared in both species)
			if human_subfamily in subfamily_toClass:
				if mouse_subfamily in subfamily_toClass:
					human_class = subfamily_toClass[human_subfamily]
					mouse_class = subfamily_toClass[mouse_subfamily]
					if human_class == mouse_class:
						orthologous_TE_numbers[1] += 1
						orthologous_TEs['human_class'].append([human_TE, human_subfamily, str(human_age)])
						orthologous_TEs['mouse_class'].append([mouse_TE, mouse_subfamily, str(mouse_age)])
						is_ortholog = True
						#See if orthologous TEs are associated with the same cCRE types
						orthologous_TE_ccres['total']['all'] += 1
						ccre_mouse = mouse_cCREs[mm10_ccre][0]
						orthologous_TE_ccres['total'][ccre_mouse] += 1
						if mm10_ccre in human_cCREs:
							if human_cCREs[mm10_ccre] != []:
								ccre_human = human_cCREs[mm10_ccre][0]
								if ccre_human == ccre_mouse:
									orthologous_TE_ccres['same']['all'] += 1
									orthologous_TE_ccres['same'][ccre_mouse] += 1
									o.write('\tTEortholog\tShared_cCRE\n')
								else:
									orthologous_TE_ccres['different']['all'] += 1
									orthologous_TE_ccres['different'][ccre_mouse] += 1
									o.write('\tTEortholog\tDifferent_cCRE\n')
									sys.stderr.write('\t'.join(di[mm10_ccre]) + '\n')
							else:
								#sys.stderr.write('\t'.join(di[mm10_ccre]) + '\n')
								orthologous_TE_ccres['mouse_only']['all'] += 1
								orthologous_TE_ccres['mouse_only'][ccre_mouse] += 1
								o.write('\tTEortholog\tMouse_cCRE\n')
						else:
							orthologous_TE_ccres['mouse_only']['all'] += 1
							orthologous_TE_ccres['mouse_only'][ccre_mouse] += 1
							o.write('\tTEortholog\tMouse_cCRE\n')
			if not is_ortholog:
				syntenic_nonOrtholog_TE['human'].append([human_subfamily, str(human_age)])
				syntenic_nonOrtholog_TE['mouse'].append([mouse_subfamily, str(mouse_age)])
				o.write('\tSyntenic')
				ccre_mouse = mouse_cCREs[mm10_ccre][0]
				if mm10_ccre in human_cCREs:
					ccre_human = human_cCREs[mm10_ccre][0]
					if ccre_mouse == ccre_human:
						o.write('\tShared_cCRE\n')
					else:
						o.write('\tDifferent_cCRE\n')
				else:
					o.write('\tMouse_cCRE\n')
				list_nonHuman_TEs.append(mm10_ccre)
		#Get non-TE background
		elif mm10_ccre not in mouse_TEs:
			if mm10_ccre in syntenic_regions_set:
				syntenic_nonTE_numbers['total']['all'] += 1
				ccre_mouse = mouse_cCREs[mm10_ccre][0]
				syntenic_nonTE_numbers['total'][ccre_mouse] += 1
				if mm10_ccre in human_cCREs:
					if human_cCREs[mm10_ccre] != []:
						ccre_human = human_cCREs[mm10_ccre][0]
						if ccre_human == ccre_mouse:
							syntenic_nonTE_numbers['same']['all'] += 1
							syntenic_nonTE_numbers['same'][ccre_mouse] += 1
							o.write('\tSyntenic\tShared_cCRE\n')
						else:
							syntenic_nonTE_numbers['different']['all'] += 1
							syntenic_nonTE_numbers['different'][ccre_mouse] += 1
							o.write('\tSyntenic\tDifferent_cCRE\n')
							#sys.stderr.write('\t'.join(di[mm10_ccre]) + '\n')
					else:
						#sys.stderr.write('\t'.join(di[mm10_ccre]) + '\n')
						syntenic_nonTE_numbers['mouse_only']['all'] += 1
						syntenic_nonTE_numbers['mouse_only'][ccre_mouse] += 1
						o.write('\tSyntenic\tMouse_cCRE\n')
				else:
					syntenic_nonTE_numbers['mouse_only']['all'] += 1
					syntenic_nonTE_numbers['mouse_only'][ccre_mouse] += 1
					o.write('\tSyntenic\tMouse_cCRE\n')
			else:
				o.write('\tMouse_specific\tMouse_cCRE\n')
		else:
			o.write('\tMouse_TE')
			if mm10_ccre in syntenic_regions_set:
				if mm10_ccre in human_cCREs:
					ccre_mouse = mouse_cCREs[mm10_ccre][0]
					ccre_human = human_cCREs[mm10_ccre][0]
					if ccre_mouse == ccre_human:
						o.write('\tShared_cCRE\n')
					else:
						o.write('\tDifferent_cCRE\n')
				else:
					o.write('\tMouse_cCRE\n')
			else:
				o.write('\tMouse_cCRE\n')

#Quantify TE contribution to species-specific (novel) cCREs
for mm10_ccre in list_mm10_ccres:
	mouse_ccre = mouse_cCREs[mm10_ccre][0]
	#Get mouse TE-cCRE numbers
	if mm10_ccre in mouse_TEs:
		if mm10_ccre not in human_cCREs: #Not cCRE in human, i.e. novel
			if mm10_ccre in list_nonHuman_TEs: #Not orthologous with a human TE
				novel_cCREs['TE']['all'] += 1
				novel_cCREs['TE'][mouse_ccre] += 1
			else:
				novel_cCREs['TE']['all'] += 1
				novel_cCREs['TE'][mouse_ccre] += 1
		else:
			shared_cCREs['TE']['all'] += 1
			with open(sys.argv[-1] + '_shared_cCRE_TE1species', 'a') as o:
				o.write('\t'.join(di[mm10_ccre]) + '\n')
	#Get mouse non-TE cCRE numbers
	else:
		if mm10_ccre not in human_cCREs: #Not cCRE in human, i.e. novel
			novel_cCREs['non-TE']['all'] += 1
			novel_cCREs['non-TE'][mouse_ccre] += 1
		else:
			shared_cCREs['non-TE']['all'] += 1

#Write orthologous TE ages to another output file
with open(sys.argv[-1] + '_orthologTE_ages', 'w') as o:
	o.write('Species\tTE_coords\tSubfamily\tSubfamily_age\n')
	for i in range(len(orthologous_TEs['human'])):
		o.write('Human\t' + '\t'.join(orthologous_TEs['human'][i]) + '\n')
		o.write('Mouse\t' + '\t'.join(orthologous_TEs['mouse'][i]) + '\n')
with open(sys.argv[-1] + '_orthologTE_ages_relatedSubfamilies', 'w') as o:
	o.write('Species\tTE_coords\tSubfamily\tSubfamily_age\n')
	for i in range(len(orthologous_TEs['human_class'])):
		o.write('Human\t' + '\t'.join(orthologous_TEs['human_class'][i]) + '\n')
		o.write('Mouse\t' + '\t'.join(orthologous_TEs['mouse_class'][i]) + '\n')

#Write syntenic, but non-orthologous TEs to output file
with open(sys.argv[-1] + '_syntenic_nonOrthologTEs', 'w') as o:
	o.write('Human_TE\tHuman_age\tMouse_TE\tMouse_age\n')
	for i in range(len(syntenic_nonOrtholog_TE['human'])):
		o.write('\t'.join(syntenic_nonOrtholog_TE['human'][i]) + '\t' + '\t'.join(syntenic_nonOrtholog_TE['mouse'][i]) + '\n')

#Write summary/statistics to standard output
print('#Category\tNumber\tPerc_all\tPerc_syntenic')
print('#Mouse_cCREs\t' + str(len(mouse_cCREs)) + '\t100\tNA')
print('#Mouse_TE-cCREs\t' + str(len(mouse_TEs)) + '\t' + str(round(len(mouse_TEs)/len(mouse_cCREs)*100,2)) + '\tNA')
print('#Human_syntenic\t' + str(len(syntenic_regions)) + '\t' + str(round(len(syntenic_regions)/len(mouse_cCREs)*100,2)) + '\t100')
print('#Human_cCREs\t' + str(len(human_cCREs)) + '\t' + str(round(len(human_cCREs)/len(mouse_cCREs)*100,2)) + '\t' + str(round(len(human_cCREs)/len(syntenic_regions)*100,2)))
print('#Human_TE\t' + str(len(human_TEs)) + '\t' + str(round(len(human_TEs)/len(mouse_cCREs)*100,2)) + '\t' + str(round(len(human_TEs)/len(syntenic_regions)*100,2)) + '\n#')
##TEs in syntenic regions that are orthologous
percent1 = round(orthologous_TE_numbers[0]/orthologous_TE_numbers[2]*100,2)
percent2 = round(orthologous_TE_numbers[1]/orthologous_TE_numbers[2]*100,2)
print('#' + str(orthologous_TE_numbers[0]) + '/' + str(orthologous_TE_numbers[2]) + ' (' + str(percent1) + '%) TEs in syntenic regions between human and mouse are the same subfamily (strict orthologs)')
print('#' + str(orthologous_TE_numbers[1]) + '/' + str(orthologous_TE_numbers[2]) + ' (' + str(percent2) + '%) TEs in syntenic regions between human and mouse are the same subfamily (subfamilies in same class)')
##Orthologous TEs mouse-human cCRE comparison
print('#\n#Orthologous TEs, mouse-human cCRE comparison')
percent_TE_same = round(orthologous_TE_ccres['same']['all']/orthologous_TE_ccres['total']['all']*100, 2)
percent_TE_different = round(orthologous_TE_ccres['different']['all']/orthologous_TE_ccres['total']['all']*100, 2)
percent_TE_mouseOnly = round(orthologous_TE_ccres['mouse_only']['all']/orthologous_TE_ccres['total']['all']*100, 2)
enrich_TE_same = round((orthologous_TE_ccres['same']['all']/orthologous_TE_ccres['total']['all'])/(syntenic_nonTE_numbers['same']['all']/syntenic_nonTE_numbers['total']['all']), 2)
enrich_TE_different = round((orthologous_TE_ccres['different']['all']/orthologous_TE_ccres['total']['all'])/(syntenic_nonTE_numbers['different']['all']/syntenic_nonTE_numbers['total']['all']), 2)
enrich_TE_mouseOnly = round((orthologous_TE_ccres['mouse_only']['all']/orthologous_TE_ccres['total']['all'])/(syntenic_nonTE_numbers['mouse_only']['all']/syntenic_nonTE_numbers['total']['all']), 2)
print('#Same cCREs:\t' + str(percent_TE_same) + '% (' + str(orthologous_TE_ccres['same']['all']) + '/' + str(orthologous_TE_ccres['total']['all']) + ')\t' + str(enrich_TE_same))
print('#Different cCREs:\t' + str(percent_TE_different) + '% (' + str(orthologous_TE_ccres['different']['all']) + '/' + str(orthologous_TE_ccres['total']['all']) + ')\t' + str(enrich_TE_different))
print('#Mouse-only cCREs:\t' + str(percent_TE_mouseOnly) + '% (' + str(orthologous_TE_ccres['mouse_only']['all']) + '/' + str(orthologous_TE_ccres['total']['all']) + ')\t' + str(enrich_TE_mouseOnly))
##Orthologous TEs mouse-human cCRE comparison - separated by cCRE type
list_ccre_types = ['PLS', 'pELS', 'dELS', 'DNase-H3K4me3', 'CTCF-only', 'DNase-only']
list_compare_types = ['same', 'different', 'mouse_only']
print('#\n#Annotation\tComparison\tcCRE_type\tEnrich_comparison\tNumber')
for compare_type in list_compare_types:
	for ccre_type in list_ccre_types:
		try:
			enrich_compare = round((orthologous_TE_ccres[compare_type][ccre_type]/orthologous_TE_ccres['total'][ccre_type])/(syntenic_nonTE_numbers[compare_type][ccre_type]/syntenic_nonTE_numbers['total'][ccre_type]), 2)
			line = '#TE\t' + compare_type + '\t' + ccre_type + '\t' + str(enrich_compare) + '\t' + str(orthologous_TE_ccres[compare_type][ccre_type])
			print(line)
		except ZeroDivisionError:
			continue
print('#\nAnnotation\tComparison\tcCRE_type\tPerc_comparison\tNumber')
for compare_type in list_compare_types:
	for ccre_type in list_ccre_types:
		try:
			percent_compare = round(orthologous_TE_ccres[compare_type][ccre_type]/orthologous_TE_ccres['total'][ccre_type]*100, 2)
			line = 'TE\t' + compare_type + '\t' + ccre_type + '\t' + str(percent_compare) + '\t' + str(orthologous_TE_ccres[compare_type][ccre_type])
			print(line)
		except ZeroDivisionError:
			continue
##Syntenic non-TE background
percent_nonTE_same = round(syntenic_nonTE_numbers['same']['all']/syntenic_nonTE_numbers['total']['all']*100, 2)
percent_nonTE_different = round(syntenic_nonTE_numbers['different']['all']/syntenic_nonTE_numbers['total']['all']*100, 2)
percent_nonTE_mouseOnly = round(syntenic_nonTE_numbers['mouse_only']['all']/syntenic_nonTE_numbers['total']['all']*100, 2)
print('#\n#Non-TE and syntenic, human-mouse cCRE comparison')
print('#Same cCREs:\t' + str(percent_nonTE_same) + '% (' + str(syntenic_nonTE_numbers['same']['all']) + '/' + str(syntenic_nonTE_numbers['total']['all']) + ')')
print('#Different cCREs:\t' + str(percent_nonTE_different) + '% (' + str(syntenic_nonTE_numbers['different']['all']) + '/' + str(syntenic_nonTE_numbers['total']['all']) + ')')
print('#Mouse-only cCREs:\t' + str(percent_nonTE_mouseOnly) + '% (' + str(syntenic_nonTE_numbers['mouse_only']['all']) + '/' + str(syntenic_nonTE_numbers['total']['all']) + ')')
##Syntenic non-TE background - separated by cCRE type
print('#\n#Annotation\tComparison\tcCRE_type\tPerc_comparison\tNumber')
for compare_type in list_compare_types:
	for ccre_type in list_ccre_types:
		try:
			percent_compare = round(syntenic_nonTE_numbers[compare_type][ccre_type]/syntenic_nonTE_numbers['total'][ccre_type]*100, 2)
			line = 'non-TE\t' + compare_type + '\t' + ccre_type + '\t' + str(percent_compare) + '\t' + str(syntenic_nonTE_numbers[compare_type][ccre_type])
			print(line)
		except ZeroDivisionError:
			continue

##Novel cCREs in mouse compared to human
print('#\n#Novel cCRE numbers')
print('#Annotation\tcCRE_type\tPerc_TE\tPerc_cCRE\tEnrich_TE/all\tNumber')
for ccre_type in list_ccre_types:
	try:
		percent_te = round(novel_cCREs['TE'][ccre_type]/novel_cCREs['TE']['all']*100, 2)
		percent_all = round(novel_cCREs['TE'][ccre_type]/(novel_cCREs['TE'][ccre_type] + novel_cCREs['non-TE'][ccre_type])*100, 2)
		enrich_te = round((novel_cCREs['TE'][ccre_type]/novel_cCREs['TE']['all'])/(novel_cCREs['non-TE'][ccre_type]/novel_cCREs['non-TE']['all']), 2)
		line = '#TE\t' + ccre_type + '\t' + str(percent_te) + '\t' + str(percent_all) + '\t' + str(enrich_te) + '\t' + str(novel_cCREs['TE'][ccre_type])
		print(line)
	except ZeroDivisionError:
		continue
print('#\n#Annotation\tcCRE_type\tPerc_nonTE\tPerc_cCRE\tNumber')
for ccre_type in list_ccre_types:
	try:
		percent_nonTE = round(novel_cCREs['non-TE'][ccre_type]/novel_cCREs['non-TE']['all']*100, 2)
		percent_all = round(novel_cCREs['non-TE'][ccre_type]/(novel_cCREs['TE'][ccre_type] + novel_cCREs['non-TE'][ccre_type])*100, 2)
		line = '#non-TE\t' + ccre_type + '\t' + str(percent_nonTE) + '\t' + str(percent_all) + '\t' + str(novel_cCREs['non-TE'][ccre_type])
		print(line)
	except ZeroDivisionError:
		continue

print('\n#All novel TE cCREs:\t' + str(novel_cCREs['TE']['all']))
print('#All shared TE cCREs:\t' + str(shared_cCREs['TE']['all']))
print('#All novel non-TE cCREs:\t' + str(novel_cCREs['non-TE']['all']))
print('#All shared non-TE cCREs:\t' + str(shared_cCREs['non-TE']['all']))

