'''
Get enriched motifs from AME using a stringent criteria.
Combines enriched motifs from AME together based on their motif archetype cluster. Assumes using HOCOMOCOv11 human database.
Picks the most significant motif for each cluster to represent.
Usage: python3 get_stringent_enriched_motifs.py <input AME tsv file> <path to archetype_motif_info.txt> <path to archetype_clusters_info.txt> <output file>
'''

import sys

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

order = []
di = {}
non_clusters = []
with open(sys.argv[1], 'r') as f:
	header = f.readline().rstrip().split('\t')
	new_header = 'Cluster_id\tCluster_name\tTop_motif\tAll_motifs\t' + '\t'.join(header[5:]) + '\n'
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
			if cluster_name not in order:
				order.append(cluster_name)
				di[cluster_name] = [cluster_id, cluster_name, motif, motif] + fields[5:]
			else:
				di[cluster_name][3] = di[cluster_name][3] + ',' + motif
		else: #Motif id not part of a motif archetype
			non_clusters.append(motif)
			di[motif] = [cluster_id, cluster_name, motif, motif] + fields[5:]

with open(sys.argv[-1], 'w') as o:
	o.write(new_header)
	for cluster_name in order:
		num_cluster_enriched = len(di[cluster_name][3].split(','))
		if (num_cluster_enriched / motifs_in_archetype[cluster_name]) >= 0.5: #Consider motif archetype as enriched if at least 50% of motifs in archetype are enriched
			o.write('\t'.join(di[cluster_name]) + '\n')
	for motif in non_clusters:
		o.write('\t'.join(di[motif]) + '\n')

