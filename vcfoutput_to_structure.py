import getopt, csv, string, sys

outputfile = './data.str'

#parse parameters

try:
	opts, args = getopt.getopt(sys.argv[1:],"hi:o:")
except getopt.GetoptError:
	print('vcfoutput_to_structure.py -i <inputfile> -o <outputfile>')
	print('More details: vcfoutput_to_structure.py -h')
	sys.exit(2)
for opt, arg in opts:
	if opt == '-h':
		print('Script for creating consensus sequences for a species or an entire alignment, across a folder containing multiple fasta files.')
		print('\nUse: vcfoutput_to_structure.py -i <inputfile> -o <outputfile>')
		print('\nParameters')
		print('h: Help, i.e. get this information')
		print('i: Input file of SNP data in the 012 format as produced by vcftools')
		print('o: Output file name, which will be in STRUCTURE format. Default is data.str in current working directory')
		sys.exit()
	elif opt in ("-i"):
		inputfile = arg
	elif opt in ("-o"):
		outputfile = arg

with open(inputfile, 'r') as f:
	data = list(list(rec) for rec in csv.reader(f, delimiter='\t'))
f.close()

# data[0] is a row = one sample; data[x][0] is the sample number, other columns are 0/1/2/-1, with 1 for heterozygous and -1 for missing; input file has no header
# output format needs to be two rows per sample with 0 for allele absent, 1 for allele present, and -9 for missing data; first six columns are ignored, rest are data

#write output file

outfile = open(outputfile, "w")
for sample in data:
	outfile.write('-9\t-9\t-9\t-9\t-9\t-9')   # six columns of nonsense
	# write row for first allele, 1 if 0 (homozygous for this allele) or 1 (heterozygous), 0 if 2 (homozygous for other allele)
	for allele in sample[1:]:
		if (allele == '0') or (allele == '1'):
			outfile.write('\t0')
		else:
			if (allele == '2'):
				outfile.write('\t1')
			else:
				outfile.write('\t-9')
	outfile.write("\n")
	outfile.write('-9\t-9\t-9\t-9\t-9\t-9')   # six columns of nonsense
	# write row for second allele, 1 if 2 (homozygous for this allele) or 1 (heterozygous), 0 if 0 (homozygous for other allele)
	for allele in sample[1:]:
		if (allele == '2') or (allele == '1'):
			outfile.write('\t1')
		else:
			if (allele == '0'):
				outfile.write('\t0')
			else:
				outfile.write('\t-9')
	outfile.write("\n")
outfile.close()
