'''
Parses output of Needle alignment of TE subfamily consensus sequences and individual elements.
Prints the alignment from Needle in fasta format to standard output.
Usage: parse_needleOutput_fasta.py <Needle output file>
'''

import sys

if len(sys.argv) != 2:
	sys.exit(__doc__)

with open(sys.argv[1], 'r') as f:
	name = ""
	seq = ""
	consensus = ""
	score = 0
	for line in f:
		if line.strip() == "":
			continue
		if line.rstrip('\n').replace(' ', '').replace('|', '').replace('.', '') == "": #if the line is only indicating base match/mismatches, skip it
			continue
		if line.startswith('#'):
			if "# 2:" in line:
				fasta_name = line.rstrip().split()[-1]
				continue
			elif "?" not in line and "chr" not in line:
				continue
		if "consensus" in line:
			fields = line.rstrip().split()
			consensus += fields[2]
			continue
		fields = line.rstrip().split()
		#Check to make sure the name of the element is the same
		if name == "": #if the first element line was reached, get the name
			name = fields[0]
		else:
			if name != fields[0]: #if not the first line, check that the name for the line is same
				print(name + '\t' + sys.argv[1])
				print(fields[0])
				print('Something went wrong with parsing the Needle alignment. More than one sequence not labeled as consensus.')
				sys.exit()
		#Add the alignment sequences on separate lines together
		seq += fields[2]

#Write alignment to standard output
print('>' + fasta_name)
print(seq)
print(consensus)

