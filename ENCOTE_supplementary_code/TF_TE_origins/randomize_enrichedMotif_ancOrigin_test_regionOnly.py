'''
Randomly selects motifs to get a background distribution for motif ancestral origin.
The number of randomly selected motifs depends on the number of enriched motif archetypes for the TE subfamily.
This version is adjusted for "region only" analysis, which only counts motifs as ancestral/consensus-derived if the motif is in the same general area in the consensus sequence.
Runs 1000 tests by default (edit script to change).
Prints the results of the randomizations to standard output. Summaries are written to <output file>.
Usage: python3 randomize_enrichedMotif_ancOrigin_test_regionOnly.py <observed enriched motif origin file> <consensus motifs file> <archetype_motif_info.txt> <motif location info directory> <output file>
'''

import sys, os
import random
import statistics as stats

num_randoms = 1000 #1000 randomization tests, change this if more or fewer randomizations are desired

if len(sys.argv) != 6:
	sys.exit(__doc__)

subfamilies = []
observed_motifs = {}
subfamily_info = {}
with open(sys.argv[1], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		subfamilies.append(subfamily)
		subfamily_info[subfamily] = fields[:4]
		num_motifs = int(fields[3])
		observed_motifs[subfamily] = num_motifs

all_motifs = []
consensus_motifs = {}
consensus_motif_starts = {}
length_motif = {}
with open(sys.argv[2], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		motif = fields[1]
		starts = [(int(x)-1) for x in fields[2].split(',')] #Convert start positions to 0-based
		if motif not in all_motifs:
			all_motifs.append(motif)
			#Get motif length
			stops = [(int(x)-1) for x in fields[3].split(',')]
			length_motif[motif] = stops[0] - starts[0] + 1
		if subfamily not in consensus_motifs:
			consensus_motifs[subfamily] = []
			consensus_motif_starts[subfamily] = {}
		consensus_motifs[subfamily].append(motif)
		#Get consensus motif start positions
		consensus_motif_starts[subfamily][motif] = starts

archetype_info = {}
with open(sys.argv[3], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		cluster_id = fields[0]
		motif = fields[1]
		archetype_info[motif] = cluster_id

motif_locations = {}
motif_locations_dir = sys.argv[4].rstrip('/') + '/'
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
motif_anc_origin_perc = {}
for subfamily in motif_locations:
	if subfamily not in consensus_motif_starts: #Subfamily has not consensus motif start positions, skip
		continue
	for motif in motif_locations[subfamily]:
		if motif not in consensus_motif_starts[subfamily]: #Subfamily consensus does not have motif, set percent ancestral origin to 0 and go to next motif
			motif_anc_origin_perc[subfamily][motif] = 0
			continue
		anc_motifs = 0
		total_motifs = 0
		for motif_start in motif_locations[subfamily][motif]:
			for consensus_start in consensus_motif_starts[subfamily][motif]:
				if ((consensus_start - 10) <= motif_start) and ((consensus_start + 10) >= motif_start): #Element motif start is found at the consensus motif +/- wiggle room
					anc_motifs += 1
					break
			total_motifs += 1
		try:
			ancestral_prop = anc_motifs/total_motifs
		except ZeroDivisionError:
			ancestral_prop = 0
		if subfamily not in motif_anc_origin_perc:
			motif_anc_origin_perc[subfamily] = {}
		motif_anc_origin_perc[subfamily][motif] = ancestral_prop

random_results = {}
for subfamily in subfamilies:
	random_results[subfamily] = []
	for random_run in range(num_randoms):
		num_motifs = observed_motifs[subfamily]
		random_motifs = []
		random_archetypes = []
		while len(random_motifs) < num_motifs:
			motif = random.choice(all_motifs)
			if motif in archetype_info:
				archetype = archetype_info[motif]
			else:
				continue
			if archetype not in random_archetypes:
				random_motifs.append(motif)
				random_archetypes.append(archetype)
		perc_random_in_consensus = 0
		for motif in random_motifs:
			if motif in motif_anc_origin_perc[subfamily]:
				perc_random_in_consensus += motif_anc_origin_perc[subfamily][motif]
		percent = perc_random_in_consensus/len(random_motifs)*100
		print(subfamily + '\t' + str(percent))
		random_results[subfamily].append(percent)

with open(sys.argv[-1], 'w') as o:
	o.write('Subfamily\tClass\tKimura_div\tNum_archetypes\tNum_randoms\tRandom_perc_anc\tMedian_random_perc_anc\n')
	for subfamily in subfamilies:
		o.write('\t'.join(subfamily_info[subfamily]) + '\t' + str(num_randoms) + '\t' + str(stats.mean(random_results[subfamily])) + '\t' + str(stats.median(random_results[subfamily])) + '\n')

