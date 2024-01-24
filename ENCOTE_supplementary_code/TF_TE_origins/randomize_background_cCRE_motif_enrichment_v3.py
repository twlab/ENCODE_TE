'''
Randomizes the selected background elements for motif enrichment.
Background elements are randomly selected while keeping length distribution equivalent to foreground elements.
	Outputs each randomization to output file based on basename (each randomization is basename + number).
This version (v3) tries to maximize the number of elements in background while keeping the length distribution the same.
	Uses KS-test to help decide when number of background elements has been maximized.
	Also uses KS-test to determine whether background length distribution can even be similar to foreground length distribution.
	Difference with v2 is using 0.01 as hard p-value cutoff for KS-test.
Usage: python3 randomize_background_cCRE_motif_enrichment_v3.py <input foreground fasta sequences> <input background fasta sequences> <number of randomizations> <output file basename>
'''

import sys
import math
import random
import statistics as stats
from scipy.stats import ks_2samp

if len(sys.argv) != 5:
	sys.exit(__doc__)

#KS-test p-value cutoff
cutoff = 0.01

try:
	num_randoms = int(sys.argv[3])
except TypeError:
	sys.stderr.write('<number of randomizations> must be an integer')
	sys.exit(__doc__)

foreground_di = {}
foreground_di['lengths'] = []
min_len = 1000000
max_len = 0
with open(sys.argv[1], 'r') as f:
	for line in f:
		if line.startswith('>'):
			name = line.rstrip('\n').lstrip('>')
			if name == 'lengths':
				sys.stderr.write("foreground sequences has name 'lengths' which conflicts with script. Changes as necessary.")
				sys.exit()
			seq = f.readline().rstrip('\n')
			length = len(seq)
			foreground_di[name] = length
			foreground_di['lengths'].append(length)
			if length < min_len:
				min_len = length
			if length > max_len:
				max_len = length
num_foreground = len(foreground_di) - 1
if (num_foreground) == 0: #No foreground elements
	sys.stderr.write('Foreground: ' + sys.argv[1] + '\t' + 'Background: ' + sys.argv[2] + '\n')
	sys.stderr.write('No foreground elements. No randomized background needed (no enrichment by default).\n')
	sys.exit()

background_di = {}
background_di['lengths'] = {}
with open(sys.argv[2], 'r') as f:
	for line in f:
		if line.startswith('>'):
			name = line.rstrip('\n').lstrip('>')
			if name == 'lengths':
				sys.stderr.write("background sequences has name 'lengths' which conflicts with script. Changes as necessary.")
				sys.exit()
			seq = f.readline().rstrip('\n')
			length = len(seq)
			if length < min_len:
				continue
			if length > max_len:
				max_len = length
			background_di[name] = seq
			background_di['lengths'][name] = length

num_bins = math.ceil((len(background_di) - 1) / 100) #Try average 100 background elements per bin
try:
	if (max_len - min_len) > 1000: #Over 1000bp difference in length between longest and shortest foreground elements
		binwidth = max([50, (max_len - min_len)/num_bins]) #Either 50bp bins or custom bin width depending on above number of bins
	else: #Shorter elements
		binwidth = max([10, (max_len - min_len)/num_bins]) #Either 10bp bins or custom bin width depending on above number of bins
except ZeroDivisionError:
	sys.stderr.write('Foreground: ' + sys.argv[1] + '\t' + 'Background: ' + sys.argv[2] + '\n')
	sys.stderr.write('No background elements. No possible randomized background.\n')
	sys.exit()
if binwidth == 50:
	num_bins = math.ceil((max_len - min_len) / 50)
elif binwidth == 10:
	num_bins = math.ceil((max_len - min_len) / 10)
bin_boundaries = {}
for i in range(num_bins):
	bin_count = 'bin' + str(i+1)
	bin_boundaries[bin_count] = min_len + binwidth*(i+1)

foreground_bins = {}
for name in foreground_di:
	if name == 'lengths':
		continue
	length = foreground_di[name]
	for i in range(num_bins):
		bin_count = 'bin' + str(i+1)
		if length <= bin_boundaries[bin_count]:
			if bin_count not in foreground_bins:
				foreground_bins[bin_count] = 0
			foreground_bins[bin_count] += 1
			break

background_binned = {}
for name in background_di['lengths']:
	length = background_di['lengths'][name]
	for i in range(num_bins):
		bin_count = 'bin' + str(i+1)
		if length <= bin_boundaries[bin_count]:
			if bin_count not in background_binned:
				background_binned[bin_count] = []
			background_binned[bin_count].append(name)
			break

def deal_with_problem_bins(bin_count, num_elements, used_elements):
	selected_elements = []
	added_num = 0
	if bin_count in background_binned: #If any background elements in this bin, add it
		selected_elements += [element for element in background_binned[bin_count] if element not in used_elements]
		added_num += len(selected_elements)
	bin_num = int(bin_count.strip('bin'))
	#Add elements to random pool from other bins
	adjust = 1
	if bin_num == num_bins: #At last bin, try previous bins
		possible_pool = []
		while len(possible_pool) < (num_elements-added_num): #Get enough elements from the previous bins for random selection (not previously selected)
			bin_adjust = 'bin' + str(bin_num-adjust)
			if bin_adjust in background_binned:
				possible_pool += [element for element in background_binned[bin_adjust] if element not in used_elements]
			adjust += 1
		selected_elements += random.sample(possible_pool, (num_elements-added_num))
	elif bin_num == 0: #At first bin, try next bins
		possible_pool = []
		while len(possible_pool) < (num_elements-added_num): #Get enough elements from the previous bins for random selection (not previously selected)
			bin_adjust = 'bin' + str(bin_num+adjust)
			if bin_adjust in background_binned:
				possible_pool += [element for element in background_binned[bin_adjust] if element not in used_elements]
			adjust += 1
		selected_elements += random.sample(possible_pool, (num_elements-added-num))
	else: #If bin is in the middle somewhere, take the bins around it
		possible_pool_high = []
		possible_pool_low = []
		possible_pool_total = 0
		while possible_pool_total < (num_elements-added_num): #Get enough elements from the surround bins for random selection (not previously selected)
			bin_adjust = 'bin' + str(bin_num+adjust)
			if bin_adjust in background_binned:
				possible_pool_high += [element for element in background_binned[bin_adjust] if element not in used_elements]
			bin_adjust = 'bin' + str(bin_num-adjust)
			if bin_adjust in background_binned:
				possible_pool_low += [element for element in background_binned[bin_adjust] if element not in used_elements]
			possible_pool_total = len(possible_pool_high) + len(possible_pool_low)
			adjust += 1
		if len(possible_pool_high) < (int((num_elements-added_num)/2)): #Higher bin has less than half needed number of elements 
			selected_elements += possible_pool_high
			selected_elements += random.sample(possible_pool_low, (num_elements - added_num - len(possible_pool_high)))
		elif len(possible_pool_low) < (int((num_elements-added_num)/2)): #Lower bin has less than half needed number of elements
			selected_elements += possible_pool_low
			selected_elements += random.sample(possible_pool_high, (num_elements - added_num - len(possible_pool_low)))
		else: #Enough higher and lower bin elements
			if len(possible_pool_low) >= len(possible_pool_high):
				try:
					selected_elements += random.sample(possible_pool_high, ((num_elements-added_num) - int((num_elements-added_num)/2)))
					selected_elements += random.sample(possible_pool_low, int((num_elements-added_num)/2))
				except ValueError:
					selected_elements += random.sample(possible_pool_high, int((num_elements-added_num)/2))
					selected_elements += random.sample(possible_pool_low, ((num_elements-added_num) - int((num_elements-added_num)/2)))
			else:
				try:
					selected_elements += random.sample(possible_pool_low, ((num_elements-added_num) - int((num_elements-added_num)/2)))
					selected_elements += random.sample(possible_pool_high, int((num_elements-added_num)/2))
				except ValueError:
					selected_elements += random.sample(possible_pool_low, int((num_elements-added_num)/2))
					selected_elements += random.sample(possible_pool_high, ((num_elements-added_num) - int((num_elements-added_num)/2)))
	return selected_elements

out_basename = sys.argv[-1]
num_random_tries = 0
random_iter = 0
while random_iter < num_randoms:
	random_elements = {}
	used_elements = {}
	if (num_foreground) > (len(background_di)-1): #More foreground elements than background elements after initial length filter, take all background elements
		for element in background_di:
			if element != 'lengths':
				used_elements[element] = 'Y'
	else:
		problem_bins = []
		for i in range(num_bins):
			bin_count = 'bin' + str(i+1)
			if bin_count in foreground_bins: #There are foreground elements in this bin
				num_elements = foreground_bins[bin_count]
				if bin_count in background_binned: #There are background elements in this bin
					try:
						random_elements[bin_count] = random.sample(background_binned[bin_count], num_elements)
					except ValueError: #More foreground elements than background elements in this bin
						problem_bins.append(bin_count)
						continue
					for element in random_elements[bin_count]:
						used_elements[element] = 'Y'
				else: #No background elements in this bin
					problem_bins.append(bin_count)
			else: #No foreground elements in this bin
				continue
		for bin_count in problem_bins:
			if bin_count not in foreground_bins:
				continue
			#Get the necessary elements for the bin by taking unused (not yet selected) elements from flanking bins
			random_elements[bin_count] = deal_with_problem_bins(bin_count, foreground_bins[bin_count], used_elements) 
			for element in random_elements[bin_count]:
				used_elements[element] = 'Y'
	#Check if randomly selected background is different from foreground
	background = []
	for element in used_elements:
		background.append(background_di['lengths'][element])
	stat, pval = ks_2samp(foreground_di['lengths'], background, alternative='less')
	if pval < cutoff:
		num_random_tries += 1
		if num_random_tries >= num_randoms and (random_iter/num_random_tries) < 0.5: #Tried random sampling at least as many times as specified and less than 50% success rate
			subfamily = sys.argv[1].split('/')[-1].replace('_overlap.fa', '')
			sys.stderr.write(subfamily + '\tRandom selection from background does not consistently yield the same length distribution as foreground\t' + 
					str(round(random_iter/num_random_tries*100)) + '%\n')
			sys.exit()
		continue
	#Try to add more backgrounds elements while maintaining length distribution
	num_elements_remaining = (len(background_di)-1) - len(used_elements)
	max_add_cycles = 100 #Maximum number of cycles to add all remaining elements
	min_elements_per_cycle = 10 #Minimum number of elements to add per cycle (on average)
	num_add_cycles = min([max_add_cycles, math.ceil(num_elements_remaining/min_elements_per_cycle)])
	added_bins = {}
	reported = False
	for i in range(num_add_cycles):
		#Add randomly selected elements while keeping the same approximate distribution of bins
		elements_to_add = []
		for bin_count in foreground_bins:
			if bin_count not in added_bins:
				added_bins[bin_count] = 0
			num_to_add = int(((i+1)/num_add_cycles)*(foreground_bins[bin_count]/num_foreground)*num_elements_remaining - added_bins[bin_count])
			if num_to_add == 0:
				continue
			try:
				bin_elements_to_add = random.sample([x for x in background_binned[bin_count] if x not in used_elements], num_to_add)
			except ValueError:
				bin_elements_to_add = deal_with_problem_bins(bin_count, num_to_add, used_elements)
			except KeyError: #No background elements in current bin
				bin_elements_to_add = deal_with_problem_bins(bin_count, num_to_add, used_elements)
			elements_to_add += bin_elements_to_add
			added_bins[bin_count] += len(bin_elements_to_add)
		#Check if adding new elements maintains length distribution (is NOT significantly different from foreground)
		for element in elements_to_add:
			background.append(background_di['lengths'][element])
		stat, pval = ks_2samp(foreground_di['lengths'], background, alternative='less')
		if pval < cutoff: #Not multiple test correcting because not obtaining totally new samples
			#Newly added elements would cause randomly selected background to be significantly different than foreground by length distribution
			#Do NOT add new elements, and output randomly selected background
			with open(out_basename + str(random_iter+1), 'w') as o:
				for element in used_elements:
					o.write('>' + element + '\n')
					o.write(background_di[element] + '\n')
			reported = True
			break
		else:
			#Add new elements and proceed to next cycle
			for element in elements_to_add:
				used_elements[element] = 'Y'
	if reported: 
		pass
	else: #Have added all possible background elements, write to output
		with open(out_basename + str(random_iter+1), 'w') as o:
			for element in used_elements:
				o.write('>' + element + '\n')
				o.write(background_di[element] + '\n')
	num_random_tries += 1
	random_iter += 1

