'''
Selects background elements based on consensus coverage.
This version (v4) uses distance of start/end positions of alignment to choose an initial set of background elements.
	Then adds additional background elements until reaching statistical difference between foreground/background.
Usage: python3 select_coverage_background_cCRE_motif_enrichment_v4.py <alignFasta file> <foreground elements bed file> <background elements bed file> <output file basename>
'''

import sys
import math
import numpy as np
from scipy.spatial import distance_matrix
from scipy.stats import ks_2samp

if len(sys.argv) != 5:
	sys.exit(__doc__)

#Constant values, change as necessary
p_thresh = 0.05 #p-value threshold for calling significance in KS 2-sample test comparing foreground vs. background coverage
stat_thresh = 0.01 #KS 2-sample test statistic threshold (representing maximum cumulative difference between foreground and background)

foreground = {}
shortest = 1000000
with open(sys.argv[2], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		element = fields[0] + ':' + fields[1] + '-' + fields[2] + '(' + fields[5] + ')'
		foreground[element] = 'Y'
		length = int(fields[2]) - int(fields[1])
		if length < shortest:
			shortest = length
try:
	subfamily = fields[3]
except NameError:
	subfamily = ''

background = {}
with open(sys.argv[3], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		if subfamily == '': #No foreground elements, skip background after getting subfamily name
			subfamily = fields[3]
			break
		if fields[3] != subfamily:
			sys.stderr.write('Inconsistent subfamily name in ' + sys.argv[2] + ' (' + subfamily + ') and ' + sys.argv[3] + ' (' + fields[3] + ')\n')
			sys.stderr.write('\tExiting...\n')
			sys.exit()
		element = fields[0] + ':' + fields[1] + '-' + fields[2] + '(' + fields[5] + ')'
		length = int(fields[2]) - int(fields[1])
		if length < shortest:
			#background[element] = 'N'
			continue
		else:
			background[element] = 'Y'

num_foreground = len(foreground)
num_background = len(background)
if num_foreground == 0: #No foreground, report and exit
	sys.stderr.write(subfamily + '\tNo foreground elements\n')
	sys.exit()
elif num_background == 0: #No background, report and exit
	sys.stderr.write(subfamily + '\tNo background elements\n')
	sys.exit()
elif num_foreground > num_background: #More foreground than background, report and exit
	sys.stderr.write(subfamily + '\tMore foreground than background elements\n')
	sys.exit()

foreground_aligns = {}
background_aligns = []
background_elements = []
di_coverage = {}
di_seq = {}
foreground_coverage = []
background_coverage = []
try:
	with open(sys.argv[1], 'r') as f:
		for line in f:
			if line.startswith('>'):
				element = line.rstrip('\n').lstrip('>')
				seq = f.readline().rstrip('\n')
				ref = f.readline().rstrip('\n')
				if element in background:
					if background[element] == 'N': #Background element is shorter than the shortest foreground element, skip
						continue
				started = False
				end = 0
				ref_position = 0
				coverage = []
				for i in range(len(ref)):
					if ref[i] != '-':
						ref_position += 1
						if seq[i] == '-': #Gap
							coverage.append(0)
							continue
						else:
							if started:
								end = ref_position
							else:
								start = ref_position
								started = True
							coverage.append(1)
				if element in foreground:
					foreground_aligns[element] = [start, end]
					if np.any(foreground_coverage):
						foreground_coverage = np.add(foreground_coverage, coverage)
					else:
						foreground_coverage = coverage
				#elif background[element] == 'Y':
				elif element in background:
					background_aligns.append([start, end])
					background_elements.append(element)
					if np.any(background_coverage):
						background_coverage = np.add(background_coverage, coverage)
					else:
						background_coverage = coverage
				di_coverage[element] = coverage
				di_seq[element] = seq.replace('-', '').upper()
except FileNotFoundError:
	sys.stderr.write(subfamily + '\tNo consensus alignment file\n')
	sys.exit()

#Get some constants
ref_length = len(ref.replace('-', ''))

#Start with finding coverages in foreground and background
foreground_proportions = np.divide(foreground_coverage, num_foreground)
foreground_forKS = []
for i in range(len(foreground_coverage)):
	foreground_forKS.extend([i]*foreground_coverage[i])
background_forKS = []
for i in range(len(background_coverage)):
	background_forKS.extend([i]*background_coverage[i])
#Test if foreground and background coverage (proportion of elements covering each consensus base) is equal
stat, pval = ks_2samp(foreground_forKS, background_forKS)
#sys.stderr.write(str(stat) + '\t' + str(pval) + '\n')
if (pval >= p_thresh) or (abs(stat) < stat_thresh) : #No significant difference between foreground and background coverage, either p-value or max cumulative difference threshold not reached
	#Report subfamily to standard output
	#print(subfamily + '\tNo difference between foreground and background coverage')
	#Write all background elements used for KS test (all elements at least as long as the shortest foreground element) to output files
	#Write background sequences to output
	with open(sys.argv[-1] + '_background.fa', 'w') as o:
		for element in sorted(background):
			o.write('>' + element + '\n' + di_seq[element] + '\n')
	#Write background elements to output
	with open(sys.argv[-1] + '_background.bed', 'w') as o:
		for element in sorted(background):
			chrom = element.split(':')[0]
			start = element.split(':')[1].split('-')[0]
			stop = element.split('-')[1].split('(')[0]
			strand = element.split('(')[1].rstrip(')')
			o.write(chrom + '\t' + start + '\t' + stop + '\t' + subfamily + '\t0\t' + strand + '\n')
	#Write background coverage to output
	background_proportions = np.divide(background_coverage, num_background)
	with open(sys.argv[-1] + '_coverage', 'w') as o:
		o.write('foreground\t' + '\t'.join([str(x) for x in foreground_coverage]) + '\n')
		o.write('background\t' + '\t'.join([str(x) for x in background_coverage]) + '\n')
		o.write('enrichment')
		num_selected = len(background)
		for i in range(len(foreground_coverage)):
			selected = background_coverage[i]
			if selected >= 5:
				enrichment = foreground_proportions[i]/(background_proportions[i])
				o.write('\t' + str(round(enrichment,4)))
			else:
				o.write('\tNA')
		o.write('\n')
	#Exit
	sys.exit()
#Sort distances for each element with element labels
max_to_keep = math.ceil(min([len(background_aligns)/2, 80000]))
sorted_distances = {}
for element1 in foreground:
	element1_aligned = foreground_aligns[element1]
	distances = distance_matrix([element1_aligned], background_aligns)
	element_distances = []
	for i in range(len(distances[0])):
		element2 = background_elements[i]
		element_distances.append([element2, distances[0][i]])
	sorted_distances[element1] = sorted(element_distances, key=lambda element_distances: element_distances[1])[:max_to_keep]
#Select background until foreground and background coverages are no longer the same (i.e. until they are statistically different)
#First add elements in groups based on distance
num_possible = math.floor(num_background/num_foreground) #Max number of closest distances to select, based on average number of possible background elements per foreground element
#sys.stderr.write(str(num_possible) + '\n')
#sys.stderr.write(str(num_foreground) + '\n')
selected_background = {}
selected_coverage = []
for i in range(ref_length):
	selected_coverage.append(0)
to_add_coverage = selected_coverage
for i in range(num_possible):
	to_add_elements = {}
	for element1 in foreground:
		for j in range(len(sorted_distances[element1])):
			element2, distance = sorted_distances[element1][j]
			if element2 in selected_background or element2 in to_add_elements:
				continue
			else:
				to_add_elements[element2] = 'Y'
				to_add_coverage = np.add(to_add_coverage, di_coverage[element2])
				break
	#to_add_proportions = np.divide(to_add_coverage, (len(selected_background) + len(to_add_elements)))
	to_add_forKS = []
	for j in range(len(to_add_coverage)):
		to_add_forKS.extend([j]*to_add_coverage[j])
	stat, pval = ks_2samp(foreground_forKS, to_add_forKS)
	#sys.stderr.write('round' + str(i) + ':\t' + str(pval) + '\t' + str(stat) + '\n')
	#sys.stderr.write(str(foreground_coverage) + '\n')
	#sys.stderr.write(str(to_add_coverage) + '\n')
	if (pval >= p_thresh) or (abs(stat) < stat_thresh): 
		#No significant difference between foreground and selected background coverage so far, add current iteration's elements and coverage to record and continue to next iteration
		selected_coverage = to_add_coverage
		for element2 in to_add_elements:
			selected_background[element2] = 'Y'
		continue
	else:
		break

if len(selected_background) == 0: #First round of picking background elements leads to significantly different coverage compared to foreground
	sys.stderr.write(subfamily + '\tNo possible background elements set that matches coverage compared to foreground\n')
	sys.exit()
sorted_distances.clear() #Clear sorted distances dictionary once it's no longer used, should free up memory?

if len(selected_background) > 10000: #More than 10000 background elements from initial round of picking background, good enough
	pass
else: #Try to add more elements for subfamilies with fewer background elements (to increase power)
	#Next add elements until foreground/background coverage are no longer equal
	num_to_add = max([1, math.ceil((num_background-len(selected_background))/20)]) #Add 1 or 5% of remaining background at a time
	for i in range(num_background - len(selected_background)):
		test_pvals = []
		selected_forKS = []
		for j in range(len(selected_coverage)):
			selected_forKS.extend([j]*selected_coverage[j])
		for element in background:
			if background[element] == 'Y':
				if element not in selected_background:
					element_coverage = di_coverage[element]
					#test_coverage = np.add(selected_coverage, element_coverage)
					#test_proportions = np.divide(test_coverage, len(selected_background)+1)
					to_add_forKS = []
					for j in range(len(element_coverage)):
						if element_coverage[j] == 1:
							to_add_forKS.append(j)
					test_forKS = selected_forKS + to_add_forKS
					test_stat, test_pval = ks_2samp(foreground_forKS, test_forKS)
					if test_pval >= p_thresh:
						test_pvals.append([element, test_pval])
		test_pvals_sorted = sorted(test_pvals, reverse=True, key=lambda test_pvals: test_pvals[1])
		if len(test_pvals) == 0: #No elements allow foreground/background coverage to stay the same when added
			break
		if num_to_add > 1:
			test_coverage = selected_coverage
			for j in range(min([num_to_add, len(test_pvals_sorted)])):
				element, test_pval = test_pvals_sorted[j]
				element_coverage = di_coverage[element]
				test_coverage = np.add(test_coverage, element_coverage)
			test_forKS = []
			for j in range(len(test_coverage)):
				test_forKS.extend([j]*test_coverage[j])
			test_stat, test_pval = ks_2samp(foreground_forKS, test_forKS)
			if test_pval < p_thresh:
				for j in range(min([num_to_add, len(test_pvals_sorted)])):
					element, test_pval = test_pvals_sorted[j]
					selected_background[element] = 'Y'
					selected_coverage = np.add(selected_coverage, di_coverage[element])
			else: #Adding the 5% remaining background leads to high KS-test p-value, add one at a time until p-value threshold is passed
				for j in range(min([num_to_add, len(test_pvals_sorted)])):
					element, test_pval = test_pvals_sorted[j]
					element_coverage = di_coverage[element]
					to_add_forKS = []
					for j in range(len(element_coverage)):
						if element_coverage[j] == 1:
							to_add_forKS.append(j)
					test_forKS = selected_forKS + to_add_forKS
					test_stat, test_pval = ks_2samp(foreground_forKS, test_forKS)
					if test_pval >= p_thresh:
						selected_background[element] = 'Y'
						selected_coverage = np.add(selected_coverage, element_coverage)
					else:
						break
				break
		else:
			element, test_pval = test_pvals_sorted[0]
			selected_background[element] = 'Y'
			selected_coverage = np.add(selected_coverage, di_coverage[element])
		#pval = test_pval

selected_forKS = []
for j in range(len(selected_coverage)):
	selected_forKS.extend([j]*selected_coverage[j])
print('Annotation\tBase')
for base in foreground_forKS:
	print('Foreground\t' + str(base))
for base in selected_forKS:
	print('Background\t' + str(base))

#Write background sequences to output
with open(sys.argv[-1] + '_background.fa', 'w') as o:
	for element in sorted(selected_background):
		o.write('>' + element + '\n' + di_seq[element] + '\n')
#Write background elements to output
with open(sys.argv[-1] + '_background.bed', 'w') as o:
	for element in sorted(selected_background):
		chrom = element.split(':')[0]
		start = element.split(':')[1].split('-')[0]
		stop = element.split('-')[1].split('(')[0]
		strand = element.split('(')[1].rstrip(')')
		o.write(chrom + '\t' + start + '\t' + stop + '\t' + subfamily + '\t0\t' + strand + '\n')
#Write background coverage to output
with open(sys.argv[-1] + '_coverage', 'w') as o:
	o.write('foreground\t' + '\t'.join([str(x) for x in foreground_coverage]) + '\n')
	o.write('background\t' + '\t'.join([str(x) for x in selected_coverage]) + '\n')
	o.write('enrichment')
	num_selected = len(selected_background)
	for i in range(len(foreground_coverage)):
		selected = selected_coverage[i]
		if selected >= 5:
			enrichment = foreground_proportions[i]/(selected/num_selected)
			o.write('\t' + str(round(enrichment,4)))
		else:
			o.write('\tNA')
	o.write('\n')

