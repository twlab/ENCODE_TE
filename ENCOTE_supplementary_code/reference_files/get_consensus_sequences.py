'''
Get TE subfamily consensus sequences based on file from Xiaoyu (for LINEs) and RepeatMasker Library fasta file.
Requires file that connects rmsk out repeat names to rmsk align repeat names used in alignment.
Prints subfamilies without a consensus sequence to standard output.
Usage: python3 get_consensus_sequences.py <list of TE subfamilies file> <rmsk out-align connect file> <Xiaoyu's TE consensus file> <RepeatLibrary fasta file> <output file>
'''

import sys
import re
from Bio.SeqIO.FastaIO import SimpleFastaParser

if len(sys.argv) != 6:
	sys.exit(__doc__)

#Percentage of alignments required for an align repeat to be considered a true consensus sequence
align_required = 0.5

list_subfamilies = []
with open(sys.argv[1], 'r') as f:
	for line in f:
		name = line.strip()
		list_subfamilies.append(name)

connect_di = {}
with open(sys.argv[2], 'r') as f:
	header = f.readline()
	for line in f:
		fields = line.rstrip('\n').split('\t')
		rmsk_out = fields[0]
		rmsk_align = fields[1]
		percent = float(fields[-1])
		if rmsk_out not in connect_di:
			connect_di[rmsk_out] = [rmsk_align, percent]
		else:
			if percent > connect_di[rmsk_out][1]: #Replace previous rmsk align repeat name if current name has higher percentage of being associated
				connect_di[rmsk_out] = [rmsk_align, percent]

xiaoyu_consensus = {}
with open(sys.argv[3], 'r') as f:
	for name, seq in SimpleFastaParser(f):
		xiaoyu_consensus[name] = seq.upper()

replib_consensus = {}
with open(sys.argv[4], 'r') as f:
	for name, seq in SimpleFastaParser(f):
		replib_consensus[name] = seq.upper()

#Get the consensus sequence for each TE subfamily and write to output file
with open(sys.argv[-1], 'w') as o:
	for subfamily in list_subfamilies:
		LINE_match = re.match('L[0-9]', subfamily)
		if LINE_match: #If subfamily is a LINE, use subfamily name as is (i.e. do not use the align name because it is only a portion of the LINE)
			if subfamily in xiaoyu_consensus:
				if subfamily in replib_consensus: #TE subfamily in both consensus files, check agreement
					if xiaoyu_consensus[subfamily] == replib_consensus[subfamily]:
						o.write('>' + subfamily + '\n' + xiaoyu_consensus[subfamily] + '\n')
					else: #consensus files disagree, take Xiaoyu's consensus for now since it creates the full LINE rather than just a fragment of it
						o.write('>' + subfamily + '\n' + xiaoyu_consensus[subfamily] + '\n')
				else:
					o.write('>' + subfamily + '\n' + xiaoyu_consensus[subfamily] + '\n')
			elif subfamily in replib_consensus:
				o.write('>' + subfamily + '\n' + replib_consensus[subfamily] + '\n')
			else: #TE subfamily not found in either consensus file
				print(subfamily + ' not in ' + sys.argv[3] + ' or ' + sys.argv[4])
		elif subfamily in connect_di:
			align_name = connect_di[subfamily][0]
			percent = connect_di[subfamily][1]
			if percent < align_required:
				print(subfamily + '(' + align_name + ') less than ' + str(round(align_required*100)) + '% of elements')
				continue
			if align_name not in xiaoyu_consensus and align_name not in replib_consensus: #TE subfamily (alignment name) not in either consensus file, try original subfamily name
				if subfamily not in xiaoyu_consensus and subfamily not in replib_consensus: #TE subfamily not in either consensus file, report to standard output and then continue
					print(subfamily + '(' + align_name + ') not in ' + sys.argv[3] + ' or ' + sys.argv[4])
					continue
				else:
					print(subfamily + '(' + align_name + ' only) not in ' + sys.argv[3] + ' or ' + sys.argv[4])
					continue
			elif align_name in replib_consensus:
				if align_name in xiaoyu_consensus: #TE subfamily in both consensus files, look for agreement
					if replib_consensus[align_name] == xiaoyu_consensus[align_name]:
						o.write('>' + subfamily + '\n' + replib_consensus[align_name] + '\n')
					else: #consensus files disagree, take Repeat Library's consensus for now
						o.write('>' + subfamily + '\n' + replib_consensus[align_name] + '\n')
				else: #TE subfamily only in RepeatLibrary, use its consensus seuqence
					o.write('>' + subfamily + '\n' + replib_consensus[align_name] + '\n')
			elif align_name in xiaoyu_consensus: #TE subfamily on in Xiaoyu's consensus sequences (probably LINEs)
				o.write('>' + subfamily + '\n' + xiaoyu_consensus[align_name] + '\n')
		else: #TE subfamily not in rmsk out-align connect file, try using TE subfamily name as is
			if subfamily not in xiaoyu_consensus and subfamily not in replib_consensus: #TE subfamily not in either consensus file, report to standard output and then continue
				print(subfamily + ' not in ' + sys.argv[2] + ' <rmsk out-align connect file>, also not in ' + sys.argv[3] + ' or ' + sys.argv[4])
				continue
			elif subfamily in xiaoyu_consensus:
				o.write('>' + subfamily + '\n' + xiaoyu_consensus[subfamily] + '\n')
			elif subfamily in replib_consensus:
				o.write('>' + subfamily + '\n' + replib_consensus[align_name] + '\n')

