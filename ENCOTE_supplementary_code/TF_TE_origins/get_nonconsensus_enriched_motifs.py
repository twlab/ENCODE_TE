'''
Gets the number of enriched motifs (archetypes) that are not present in the subfamily consensus sequences.
Usage: python3 get_nonconsensus_enriched_motifs.py <regionOnly perMotif file> <Kimura divergence tsv file> <output file>
'''

import sys

if len(sys.argv) != 4:
	sys.exit(__doc__)

numbers = {}
order = []
with open(sys.argv[1], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		TEclass = fields[1]
		if subfamily not in numbers:
			numbers[subfamily] = [TEclass, 0, 0]
			order.append(subfamily)
		numbers[subfamily][1] += 1
		if fields[3] == '0': #Top motif of the motif archetype not found in consensus
			numbers[subfamily][2] += 1

divergence = {}
with open(sys.argv[2], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		subfamily = fields[0]
		div = fields[2]
		divergence[subfamily] = div	

with open(sys.argv[-1], 'w') as o:
	o.write('Subfamily\tClass\tKimura_div\tNum_archetypes\tNum_nonConsensus_archetypes\tPerc_nonConsensus_archetypes\n')
	for subfamily in order:
		percent = numbers[subfamily][2]/numbers[subfamily][1] * 100
		if subfamily not in divergence:
			o.write(subfamily + '\t' + numbers[subfamily][0] + '\tNA\t' + 
				str(numbers[subfamily][1]) + '\t' + str(numbers[subfamily][2]) + '\t' + str(percent) + '\n')
		else:
			o.write(subfamily + '\t' + numbers[subfamily][0] + '\t' + divergence[subfamily] + '\t' + 
				str(numbers[subfamily][1]) + '\t' + str(numbers[subfamily][2]) + '\t' + str(percent) + '\n')

