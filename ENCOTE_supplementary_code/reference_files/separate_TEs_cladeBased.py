'''
Separates TEs based on which clades they are found in.
Usage: python3 separate_TEs_cladeBased.py <list TEs file> <file to connect dfam names rmsk names> <clades to include file> <ordered clades file> <TE clades file> <output file>
'''

import sys

if len(sys.argv) != 7:
	sys.exit(__doc__)

list_TEs = []
with open(sys.argv[1], 'r') as f:
	for line in f:
		subfamily = line.rstrip()
		if subfamily not in list_TEs:
			list_TEs.append(subfamily)

connect_TEs = {}
di_dfam_names = {}
with open(sys.argv[2], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		rmsk_name = fields[0]
		dfam_name = fields[1]
		if dfam_name == 'NA':
			continue
		if rmsk_name in list_TEs:
			connect_TEs[rmsk_name] = dfam_name
			if dfam_name not in di_dfam_names:
				di_dfam_names[dfam_name] = []
			di_dfam_names[dfam_name].append(rmsk_name)

included_clades = []
with open(sys.argv[3], 'r') as f:
	for line in f:
		clade = line.rstrip()
		if clade not in included_clades:
			included_clades.append(clade)

ordered_clades = []
with open(sys.argv[4], 'r') as f:
	for line in f:
		ordered_clades.append(line.rstrip('\n').split('\t'))

di = {}
with open(sys.argv[5], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		dfam_name = fields[0]
		if dfam_name not in di_dfam_names: #No corresponding rmsk name and/or rmsk name not in the list of desired TEs
			continue
		clades = fields[1:]
		include = False
		for clade in clades:
			if clade in included_clades:
				include = True
				break
		if include:
			if len(clades) > 1:
				clade_tiers = []
				for clade in clades:
					found = False
					for i in range(len(ordered_clades)):
						for tier_clade in ordered_clades[i]:
							if clade == tier_clade:
								clade_tiers.append(i)
								found = True
								break
						if found:
							break
					if not found:
						clade_tiers.append(len(ordered_clades) + 1)
				lowest_tier = min(clade_tiers)
				if lowest_tier == (len(ordered_clades) + 1): #Clade for the subfamily is not found in the ordered clades file, stop and report
					print(dfam_name + ' with clades (' + ','.join(clades) + ') not found in <ordered clades file> ' + sys.argv[4], file=sys.stderr)
					sys.exit(__doc__)
				clade = clades[clade_tiers.index(lowest_tier)]
			else:
				clade = clades[0]
			#Check that the clade is still in included clades list (may not be if there are multiple clades listed for a subfamily)
			if clade not in included_clades:
				continue
			#Add subfamily and clade info to di
			if clade not in di:
				di[clade] = ''
			for subfamily in di_dfam_names[dfam_name]:
				di[clade] += subfamily + '\t' + clade + '\n'

with open(sys.argv[-1], 'w') as o:
	o.write('Subfamily\tClade\n')
	for clade in included_clades:
		if clade in di:
			o.write(di[clade].rstrip('\n') + '\n')

