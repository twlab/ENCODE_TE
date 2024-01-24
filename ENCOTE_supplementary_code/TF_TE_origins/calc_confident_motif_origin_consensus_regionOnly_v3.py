'''
Calculates the proportion of enriched motifs that have an ancestral origin based on presence in the TE subfamily consensus.
Only considers the region of the consensus that is aligned each element's motif (so motif presence in the consensus depends on the consensus region instead of full length).
Uses all motifs found in the subfamily (even if there were multiple motifs in the same element).
Print each TE subfamily-enriched motif archetype result to standard output.
This version uses the high confidence enriched motifs format (slightly changes which columns information is in).
Version 3 uses the top enriched motif per motif archetype that is found the consensus sequence (if any). For length and coverage controlled motifs.
Usage: python3 calc_confident_motif_origin_consensus_regionOnly_v2.py <enriched motifs directory> <consensus motifs file> <list TE subfamilies with consensus file> <TE subfamily class file> <TE subfamily Kimura distance file> <motif location info directory> <output file>
'''

import sys, os

if len(sys.argv) != 8:
	sys.exit(__doc__)

#Distance (bp) allowed for motif to still be considered ancestral/consensus derived
allowed_dist = 10 #10bp wiggle room around the consensus motif allowed

consensus_motifs = {}
with open(sys.argv[2], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		motif = fields[1]
		starts = fields[2].split(',')
		if subfamily not in consensus_motifs:
			consensus_motifs[subfamily] = {}
		consensus_motifs[subfamily][motif] = [(int(x)-1) for x in starts] #convert from 1-based to 0-based

with open(sys.argv[3], 'r') as f:
	for line in f:
		subfamily = line.strip()
		if subfamily not in consensus_motifs:
			consensus_motifs[subfamily] = []

list_classes = ["DNA", "LTR", "SINE", "LINE", "SVA", "Other"]
di_class = {}
with open(sys.argv[4], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		if "DNA" in fields[1]:
			TEclass = "DNA"
		elif "LTR" in fields[1]:
			TEclass = "LTR"
		elif "SINE" in fields[1]:
			TEclass = "SINE"
		elif "LINE" in fields[1]:
			TEclass = "LINE"
		elif "SVA" in fields[1]:
			TEclass = "SVA"
		else:
			TEclass = "Other"
		di_class[subfamily] = TEclass

divergence = {}
with open(sys.argv[5], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		div = fields[2]
		divergence[subfamily] = div

motif_locations = {}
motif_locations_dir = sys.argv[6].rstrip('/') + '/'
list_motif_location_files = os.listdir(motif_locations_dir)
for filename in list_motif_location_files:
	subfamily = ''
	motif = ''
	with open(motif_locations_dir + filename, 'r') as f:
		header = f.readline()
		for line in f:
			fields = line.rstrip('\n').split('\t')
			if fields[0] == 'Subfamily': #2nd line is actually header, this happens if no common start is found
				continue
			if subfamily == '':
				subfamily = fields[0]
			elif subfamily != fields[0]:
				sys.stderr.write(motif_locations_dir + filename + ' has ' + fields[0] + ' subfamily\n')
				sys.exit()
			if motif == '':
				motif = fields[1]
			elif motif != fields[1]:
				sys.stderr.write(motif_locations_dir + filename + ' has ' + fields[1] + ' motif\n')
				sys.exit()
			if subfamily not in motif_locations:
				motif_locations[subfamily] = {}
			if motif not in motif_locations[subfamily]:
				motif_locations[subfamily][motif] = []
			try:
				if fields[4] == 'NA':
					continue #Skip if no motif start within the consensus sequence found (happens when the element is longer than the consensus)
				else:
					motif_start = int(fields[4])
			except ValueError:
				sys.stderr.write(motif_locations_dir + filename + '\n')
				raise
			motif_locations[subfamily][motif].append(motif_start)

enriched_dir = sys.argv[1].rstrip('/') + '/'
list_files = os.listdir(enriched_dir)
subfamilies = []
di = {}
for filename in list_files:
	subfamily = filename.split('.')[0]
	subfamilies.append(subfamily)
	di[subfamily] = []
	if subfamily not in consensus_motifs:
		sys.stderr.write(subfamily + ' not in ' + sys.argv[2] + '\n')
		continue
	if subfamily not in motif_locations:
		sys.stderr.write(subfamily + ' not in ' + motif_locations_dir + '\n')
		continue
	with open(enriched_dir + filename, 'r') as f:
		header = f.readline()
		for line in f:
			fields = line.rstrip('\n').split('\t')
			cluster_id = fields[11]
			cluster_name = fields[12]
			top_motif = fields[1]
			enriched_motifs = fields[13].split(',')
			archetype_in_consensus = False
			next_motif = ''
			for motif in enriched_motifs:
				if motif in consensus_motifs[subfamily]:
					archetype_in_consensus = True
					next_motif = motif
					break
			#Used to check which top motifs are not found in motif locations directory, comment out/remove section to run script normally
			#if top_motif not in motif_locations[subfamily]:
			#	sys.stderr.write(subfamily + '_' + top_motif + ' not in ' + motif_locations_dir + '\n')
			#	continue
			if top_motif in consensus_motifs[subfamily]: #Top motif of the archetype is found in the consensus
				anc_motifs = 0
				total_motifs = 0
				if top_motif in motif_locations[subfamily]: #Top motif has a common start position relative to consensus in any element
					for motif_start in motif_locations[subfamily][top_motif]:
						for consensus_start in consensus_motifs[subfamily][top_motif]:
							if ((consensus_start - 10) <= motif_start) and ((consensus_start + 10) >= motif_start): #Element motif start is found at the consensus motif +/- wiggle room
								anc_motifs += 1
								break
						total_motifs += 1
				else: #No element has the top motif => motif is assigned 0% ancestral derived (since no motifs are found)
					pass
				try:
					ancestral_prop = anc_motifs/total_motifs
				except ZeroDivisionError:
					ancestral_prop = 0
			elif archetype_in_consensus and next_motif in consensus_motifs[subfamily]:
				anc_motifs = 0
				total_motifs = 0
				if next_motif in motif_locations[subfamily]: #Top motif has a common start position relative to consensus in any element
					for motif_start in motif_locations[subfamily][next_motif]:
						for consensus_start in consensus_motifs[subfamily][next_motif]:
							if ((consensus_start - 10) <= motif_start) and ((consensus_start + 10) >= motif_start): #Element motif start is found at the consensus motif +/- wiggle room
								anc_motifs += 1
								break
						total_motifs += 1
				else: #No element has the top motif => motif is assigned 0% ancestral derived (since no motifs are found)
					pass
				try:
					ancestral_prop = anc_motifs/total_motifs
				except ZeroDivisionError:
					ancestral_prop = 0
			elif archetype_in_consensus: #Should never get here
				sys.stderr.write(subfamily + ' ' + next_motif + ' motif has problem\n')
				sys.exit()
			else:
				ancestral_prop = 0
			if subfamily in di_class:
				TEclass = di_class[subfamily]
			else:
				TEclass = "Other"
			print(subfamily + '\t' + TEclass + '\t' + cluster_name + '\t' + str(ancestral_prop))
			di[subfamily].append([cluster_name, ancestral_prop])

with open(sys.argv[-1], 'w') as o:
	o.write('Subfamily\tClass\tKimura_div\tNum_archetypes\tPerc_anc_origin\n')
	for subfamily in subfamilies:
		if di[subfamily] == []: #No enriched motifs or no consensus, skip
			continue
		total_prop = 0
		for motif_cluster in di[subfamily]:
			total_prop += motif_cluster[1]
		total = len(di[subfamily])
		percent = total_prop/total*100
		if subfamily in di_class:
			TEclass = di_class[subfamily]
		else:
			TEclass = "Other"
		#Get the kimura distance
		if subfamily in divergence:
			div = divergence[subfamily]
		else:
			print(subfamily + ' does not have Kimura distance', file=sys.stderr)
			div = 'NA'
		o.write(subfamily + '\t' + TEclass + '\t' + div + '\t' + str(total) + '\t' + str(percent) + '\n')

