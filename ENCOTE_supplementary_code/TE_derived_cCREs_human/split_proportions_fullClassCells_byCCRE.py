'''
Splits TE-cCRE proportions by cCRE type after combining info across full classification cell types.
Usage: python3 split_proportions_fullClassCells_byCCRE.py <input dir> <output file basename>
'''

import sys, os

if len(sys.argv) != 3:
	sys.exit(__doc__)

input_dir = sys.argv[1].rstrip('/') + '/'
list_filenames = os.listdir(input_dir)

list_celltypes = []
list_TEclasses = []
list_ccre_types = []
di = {}
for filename in list_filenames:
	with open(input_dir + filename, 'r') as f:
		celltype = '_'.join(filename.split('_')[:-2])
		if celltype not in list_celltypes:
			list_celltypes.append(celltype)
		for line in f:
			if line.startswith('#'):
				continue
			fields = line.rstrip('\n').split('\t')
			if fields[0] == 'Class': #Header
				continue
			TEclass = fields[0]
			if TEclass not in list_TEclasses:
				list_TEclasses.append(TEclass)
			ccre_type = fields[1]
			if ccre_type not in list_ccre_types:
				list_ccre_types.append(ccre_type)
			if ccre_type not in di:
				di[ccre_type] = {}
			if celltype not in di[ccre_type]:
				di[ccre_type][celltype] = {}
			di[ccre_type][celltype][TEclass] = fields[2:]

output_basename = sys.argv[-1]
for ccre_type in list_ccre_types:
	filename = output_basename + '_' + ccre_type.replace(' ', '_')
	with open(filename, 'w') as o:
		o.write('cCRE_type\tCell_type\tClass\tcCRE_num\tcCRE_bp\tcCRE_num_perc\tcCRE_bp_perc\n')
		for celltype in list_celltypes:
			if celltype in di[ccre_type]:
				for TEclass in list_TEclasses:
					if TEclass in di[ccre_type][celltype]:
						o.write(ccre_type + '\t' + celltype + '\t' + TEclass + '\t' + '\t'.join(di[ccre_type][celltype][TEclass]) + '\n')

