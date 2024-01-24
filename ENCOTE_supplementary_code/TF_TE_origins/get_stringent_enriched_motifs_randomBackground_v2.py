'''
Get enriched motifs from AME using a stringent criteria.
Combines enriched motifs from AME together based on their motif archetype cluster. Assumes using HOCOMOCOv11 human database.
Picks the most significant motif for each cluster to represent.
This version consolidates AME enriched motifs across multiple randomized backgrounds.
	Version2 removes requirement for at least 50% motif archetype to be enriched.
Usage: python3 get_stringent_enriched_motifs_randomBackground_v2.py <path to randomized AME output directories> <path to archetype_motif_info.txt> <path to archetype_clusters_info.txt> <output file>
'''

import sys, os
from collections import Counter

if len(sys.argv) != 5:
	sys.exit(__doc__)

di_cluster_info = {}
with open(sys.argv[3], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		cluster_id = fields[0]
		cluster_name = fields[1]
		di_cluster_info[cluster_id] = cluster_name

di_archetypes = {}
motifs_in_archetype = {}
with open(sys.argv[2], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		cluster_id = fields[0]
		cluster_name = di_cluster_info[cluster_id]
		database = fields[2]
		if database != 'HOCOMOCO_v11': #Skip if not HOCOMOCO motif
			continue
		motif = fields[1]
		if 'HUMAN' not in motif.upper(): #Skip if not human motif (probably mouse then)
			continue
		di_archetypes[motif] = [cluster_name, cluster_id]
		if cluster_name not in motifs_in_archetype:
			motifs_in_archetype[cluster_name] = 0
		motifs_in_archetype[cluster_name] += 1

filenames = []
output_dir = sys.argv[1].rstrip('/') + '/'
for directory in os.scandir(output_dir):
	if directory.is_dir():
		if 'ame.tsv' in (os.listdir(directory.path)):
			filenames.append(directory.path + '/ame.tsv')

di = {}
cluster_counts = {}
clusters = []
non_clusters = []
for filename in filenames:
	with open(filename, 'r') as f:
		di[filename] = {}
		header = f.readline().rstrip().split('\t')
		for line in f:
			if line.startswith('#'):
				continue
			elif line.strip() == '':
				continue
			fields = line.rstrip('\n').split('\t')
			if fields[2] in di_archetypes: #First try motif_ID
				motif = fields[2]
				cluster_name = di_archetypes[motif][0]
				cluster_id = di_archetypes[motif][1]
			elif fields[3] in di_archetypes: #Next try motif_alt_ID
				motif = fields[3]
				cluster_name = di_archetypes[motif][0]
				cluster_id = di_archetypes[motif][1]
			else: #Motif id not found in archetypes
				cluster_name = 'NA'
				cluster_id = 'NA'
			if cluster_name != 'NA': #Motif id is part of a motif archetype
				if cluster_name not in clusters:
					clusters.append(cluster_name)
					cluster_counts[cluster_name] = 0
				if cluster_name not in di[filename]:
					di[filename][cluster_name] = [cluster_id, cluster_name, motif, motif]
					cluster_counts[cluster_name] += 1
				else:
					di[filename][cluster_name][3] = di[filename][cluster_name][3] + ',' + motif
			else: #Motif id not part of a motif archetype
				if motif not in non_clusters:
					non_clusters.append(motif)
					cluster_counts[motif] = 0
				di[filename][motif] = [cluster_id, cluster_name, motif, motif]
				cluster_counts[motif] += 1

with open(sys.argv[-1], 'w') as o:
	new_header = 'Cluster_id\tCluster_name\tMost_common_top_motif\tAll_top_motifs\tAll_motifs\n'
	o.write(new_header)
	for cluster_name in clusters:
		if cluster_counts[cluster_name] == len(filenames): #If all randomizations have found the motif archetype to be enriched
			top_motifs = []
			all_motifs = {}
			for filename in di:
				cluster_motifs = di[filename][cluster_name][3].split(',')
				for motif in cluster_motifs:
					if motif not in all_motifs:
						all_motifs[motif] = 0
					all_motifs[motif] += 1
				top_motifs.append(di[filename][cluster_name][2])
			most_common = Counter(top_motifs).most_common(1)[0][0]
			archetype_pass = False
			for motif in all_motifs:
				if all_motifs[motif] == len(filenames): #At least one motif is found to be enriched in all randomizations
					archetype_pass = True
					break
			if archetype_pass:
				o.write('\t'.join(di[filenames[0]][cluster_name][:2]) + '\t' + most_common + '\t' + ','.join(top_motifs) + '\t' + ','.join(all_motifs.keys()) + '\n')
		else:
			continue
	for motif in non_clusters:
		if cluster_counts[motif] == len(filenames): #If all randomizations have found the motif to be enriched
			o.write('\t'.join(di[filenames[0]][motif][:3]) + '\n')

