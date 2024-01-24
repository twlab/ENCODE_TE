'''
Uses dbSNP rs ID to find or verify genomic coordinates.
Writes GWAS entries with no easily discernible coordinates to standard output.
Writes mismatching GWAS/dbSNP coordinates to standard error.
Usage: python3 match_coords_snpID.py <input GWAS file> <dbSNP info file> <chrom sizes file> <output file>
'''

import sys

if len(sys.argv) != 5:
	sys.exit(__doc__)

sizes_di = {}
with open(sys.argv[3], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		chrom = fields[0]
		size = int(fields[1])
		try:
			if int(chrom.replace('chr', '')) <= 22: #Up to chromosome 22 (shouldn't be chromosome 0 or negative chromosomes...)
				sizes_di[chrom] = size
				continue
		except ValueError:
			pass
		if chrom == 'chrX' or chrom == 'chrY':
			sizes_di[chrom] = size
		else: #Skip non chr1-22 or chrX or chrY
			continue

snp_di = {}
with open(sys.argv[2], 'r') as f:
	for line in f:
		fields = line.rstrip('\n').split('\t')
		chrom = fields[0]
		start = int(fields[1])
		stop = int(fields[2])
		rs_id = fields[3]
		if chrom in sizes_di:
			if rs_id not in snp_di:
				snp_di[rs_id] = [[chrom, start, stop]]
			else:
				snp_di[rs_id].append([chrom, start, stop])
		else: #Skip non chr1-22 or chrX or chrY
			continue

gwas_di = {}
with open(sys.argv[1], 'r') as f:
	header = f.readline().rstrip('\n').split('\t')
	for line in f:
		fields = line.rstrip('\n').split('\t')
		pubmed_id = fields[0]
		disease_trait = fields[1]
		chrom = 'chr' + fields[2]
		if chrom != 'chr':
			if chrom in sizes_di:
				pass
			else:
				print(line.rstrip('\n'))
				continue
		position = fields[3]
		mapped_trait = fields[7]
		trait_uri = fields[8]
		if ('x' in chrom and 'x' in fields[4]) or (fields[4].count('rs') > 1 and 'x' in fields[4]):
			sys.stderr.write('SNP interaction:\t' + line)
			continue
		elif len(fields[4].split(';')) > 1:
			sys.stderr.write('Too many SNPs:\t' + line)
			continue
		elif len(fields[4].split(',')) > 1:
			sys.stderr.write('Too many SNPs:\t' + line)
			continue
		if chrom == 'chr' or position == '': #Either chromosome or position missing, should have SNP rs ID or chr:position in SNP rs ID field (in most cases...)
			if 'rs' in fields[4]:
				if fields[5] == '1':
					rs_id = 'rs' + fields[6].replace('rs', '')
				else:
					rs_id = fields[4].split('-')[0]
				if rs_id in snp_di:
					if len(snp_di[rs_id]) > 1:
						sys.stderr.write('Multiple possible coordinates:\t' + line)
						continue
					else:
						chrom, start, stop = snp_di[rs_id][0]
						continue
				else:
					sys.stderr.write('rsID not in dbSNP:\t' + line)
					continue
			elif ':' in fields[4]: #Has chr:position info (more common format)
				chrom = fields[4].replace('chr', '').replace('Chr:', '').replace('Chr', '').split(':')[0]
				if 'chr' not in chrom:
					chrom = 'chr' + chrom
				if chrom in sizes_di:
					pass
				else:
					print(line.rstrip('\n'))
					continue
				try:
					if 'chr' in fields[4].split(':')[1]:
						stop = int(fields[4].replace('chr', '').replace('Chr:', '').replace('Chr', '').replace('"', '').split(':')[2].split('-')[0]) - 1
					else:
						stop = int(fields[4].replace('chr', '').replace('Chr:', '').replace('Chr', '').replace('"', '').split(':')[1].split('-')[0]) - 1
				except ValueError:
					stop = int(fields[4].replace('Chr:', '').replace('Chr', '').replace('"', '').split(':')[1].split('-')[0].split('_')[0]) - 1
				start = stop
			elif ('chr' in fields[4] and '_' in fields[4]): #Has chr_position (uncommon format)
				chrom = fields[4].split('_')[0]
				if 'chr' not in chrom:
					chrom = 'chr' + chrom
				if chrom in sizes_di:
					pass
				else:
					print(line.rstrip('\n'))
					continue
				try:
					stop = int(fields[4].replace('chr', '').replace('Chr:', '').replace('Chr', '').split('_')[1].split('-')[0]) - 1
					start = stop
				except IndexError: #Only has "chr#" with no position, report and skip
					print(line.rstrip('\n'))
					continue
			else: #No chr:position info either, report to standard output and skip
				print(line.rstrip('\n'))
				continue
		else:
			stop = int(position) - 1
			start = stop
			#Check that start and stop of listed GWAS coords match dbSNP coords (at least fall within bounds)
			if fields[5] == '1':
				rs_id = 'rs' + fields[6].replace('rs', '')
			else:
				rs_id = fields[4].split('-')[0]
			if rs_id in snp_di:
				match = False
				match_index = 0
				for i in range(len(snp_di[rs_id])):
					snp_chrom = snp_di[rs_id][i][0]
					if chrom != snp_chrom:
						continue
					snp_start, snp_stop = snp_di[rs_id][i][1:]
					if (start >= snp_start and start <= snp_stop) and (stop >= snp_start and stop <= snp_stop):
						match = True
						match_index = i
				if not match: #No match between GWAS and dbSNP coords, report and skip
					sys.stderr.write('No match:\t' + line)
					continue
			else: #dbSNP doesn't have SNP for whatever reason, skip the coord double check step
				pass
		#Check that start and stop coords are possible given chromosome sizes
		if start > sizes_di[chrom] or stop > sizes_di[chrom]:
			sys.stderr.write('Start/stop position greater than chromosome size:\t' + line)
			continue
		#Record GWAS info
		coord = chrom + ':' + str(start) + '-' + str(stop)
		if coord not in gwas_di:
			gwas_di[coord] = {}
			gwas_di[coord][pubmed_id] = [chrom, str(start), str(stop), disease_trait, pubmed_id, mapped_trait, trait_uri]
		else:
			if pubmed_id in gwas_di[coord]: #Same variant coord, same study => duplicate variant
				continue
			else:
				gwas_di[coord][pubmed_id] = [chrom, str(start), str(stop), disease_trait, pubmed_id, mapped_trait, trait_uri]

with open(sys.argv[-1], 'w') as o:
	for coord in sorted(gwas_di):
		for pubmed_id in sorted(gwas_di[coord]):
			o.write('\t'.join(gwas_di[coord][pubmed_id]) + '\n')

