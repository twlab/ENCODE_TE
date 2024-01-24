"""-----------------------------------------------------------------------------
  * Script Name: controlMerge.v3.py
  * Description: This script will iterate through the RepeatMasker output and then
	merge together annotations that are within a distance threshold based on the 
	repLeft value of the current element. It will *only* merge elements if they 
	match up in terms of repStart, repEnd and repLeft rather than just same name
  * Created By:  Jason David Chobirko
  * Date:        April 14th, 2022

  * Thursday March 23rd, 2023: You added a check to make sure the chromosome is also
	the same before merging, since running this for the strains makes that 
	possibility happen frequently enough that it caused errors. Better late than 
	never eh?
-----------------------------------------------------------------------------"""
# First import all of the necessary modules you need
import sys, os, re

# The expected call to run the script is:
# python controlMerge.py genome.fa.out output.txt

# This opens file 1) - the output of RepeatMasker on mm10 (or any genome, really!)
file = open(sys.argv[1], 'r')

# This creates an output file for the RepeatMasker output
newFile = open(sys.argv[2], 'w')

# This takes the input file and joins all the whitespaces into a single whitespace, yay!
# The entire file object is saved here so you can iterate through this to get what you want.
test = [' '.join(x.split()) for x in file]

# Time to run the main merge method over the genome file 
# You also added a column of the divergence value so you can estimate element age
count = 1; mergeCount = 1; success = False
for x in test[3:]:
	split = re.split('[\s\n]', x)
	if (count == 1):
		preChr = split[4]; preStart = int(split[5]); preEnd = int(split[6]); preName = split[9]; preClass = split[10]; preDiv = float(split[1]); preConEnd = int(split[12])
		if (split[8] == "+"):
			preStrand = "+"; preAlign = split[11] + ":" + split[12]; preLeft = int(re.sub('[()]', '', split[13])); preConStart = int(split[11])
		else:
			preStrand = "-"; preAlign = split[12] + ":" + split[13]; preLeft = int(re.sub('[()]', '', split[11])); preConStart = int(split[13])
	else:
		curChr = split[4]; curStart = int(split[5]); curEnd = int(split[6]); curName = split[9]; curClass = split[10]; curDiv = float(split[1]); curConEnd = int(split[12])
		if (split[8] == "+"):
			curStrand = "+"; curAlign = split[11] + ":" + split[12]; curLeft = int(re.sub('[()]', '', split[13])); curConStart = int(split[11])
		else:
			curStrand = "-"; curAlign = split[12] + ":" + split[13]; curLeft = int(re.sub('[()]', '', split[11])); curConStart = int(split[13])
		if (preStrand == curStrand and preName == curName and preChr == curChr):
			# Now test if the previous and current element are reasonably adjacent to merge together via mode and left
			# This *should NOT* matter for strand since left is already calculated above. Although for + strand it should 
			# calculate left from the pre element, not cur? Old code: mode_dict[preName] - 
			if (curStrand == "+"):
				if ( (curStart - preEnd) <= (((curEnd - curStart) + (preEnd - preStart)) / 10) and ((curConStart - preConEnd) > 0) ):
				# if (((curStart - preEnd) <= (curConStart - preConEnd)) or (((curStart - preEnd) == 1) and ((curConStart - preConEnd) > 0))):
					success = True
				else: 
					success = False
			else:
				if ( (curStart - preEnd) <= (((curEnd - curStart) + (preEnd - preStart)) / 10) and ((preConStart - curConEnd) > 0) ):
				# if (((curStart - preEnd) <= (preConStart - curConEnd)) or (((curStart - preEnd) == 1) and ((preConStart - curConEnd) > 0))):
					success = True
				else: 
					success = False
			if (success):
				if (mergeCount == 1):
					tmpStart = preStart; tmpEnd = curEnd; tmpAlign = preAlign + "," + curAlign; tmpDiv = (preDiv * (preEnd - preStart)) + (curDiv * (curEnd - curStart))
				else:
					tmpEnd = curEnd; tmpAlign = tmpAlign + "," + curAlign; tmpDiv = tmpDiv + (curDiv * (curEnd - curStart))
				mergeCount += 1
			# Otherwise the strand and element are the same but they don't reasonably merge together. Output as below
			else:
				if (mergeCount > 1):
					a = newFile.write(str(preChr) + '\t' + str(tmpStart) + '\t' + str(tmpEnd) + '\t' + str(preName) + '\t' + str(preClass) + '\t' + str(preStrand) + '\t' + str(tmpAlign) + '\t' + str(round(tmpDiv / (tmpEnd - tmpStart), 2)) + '\n')
				else:
					a = newFile.write(str(preChr) + '\t' + str(preStart) + '\t' + str(preEnd) + '\t' + str(preName) + '\t' + str(preClass) + '\t' + str(preStrand) + '\t' + str(preAlign) + '\t' + str(preDiv) + '\n')
				mergeCount = 1
		else:
			if (mergeCount > 1):
				a = newFile.write(str(preChr) + '\t' + str(tmpStart) + '\t' + str(tmpEnd) + '\t' + str(preName) + '\t' + str(preClass) + '\t' + str(preStrand) + '\t' + str(tmpAlign) + '\t' + str(round(tmpDiv / (tmpEnd - tmpStart), 2)) + '\n')
			else:
				a = newFile.write(str(preChr) + '\t' + str(preStart) + '\t' + str(preEnd) + '\t' + str(preName) + '\t' + str(preClass) + '\t' + str(preStrand) + '\t' + str(preAlign) + '\t' + str(preDiv) + '\n')
			mergeCount = 1
		preChr = curChr; preStart = curStart; preEnd = curEnd; preDiv = curDiv; preName = curName; preClass = curClass; preStrand = curStrand; preAlign = curAlign; preLeft = curLeft; preConEnd = curConEnd; preConStart = curConStart
	count += 1
	if (count > (len(test) - 3)):
		if (mergeCount > 1):
			a = newFile.write(str(preChr) + '\t' + str(tmpStart) + '\t' + str(tmpEnd) + '\t' + str(preName) + '\t' + str(preClass) + '\t' + str(preStrand) + '\t' + str(tmpAlign) + '\t' + str(round(tmpDiv / (tmpEnd - tmpStart), 2)) + '\n')
		else:
			a = newFile.write(str(preChr) + '\t' + str(preStart) + '\t' + str(preEnd) + '\t' + str(preName) + '\t' + str(preClass) + '\t' + str(preStrand) + '\t' + str(preAlign) + '\t' + str(preDiv) + '\n')
newFile.close()



