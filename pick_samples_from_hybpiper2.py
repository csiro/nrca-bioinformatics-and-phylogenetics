import sys, getopt
import glob
from Bio import AlignIO
from Bio import SeqIO

#list of variables to be populated from parameters
inputfiledir = '' #-i
outputfiledir = '' #-o
samplelistfile = '' # -s
genelistfile = '' # -g

#parse parameters

try:
	opts, args = getopt.getopt(sys.argv[1:],"hi:o:s:g:")
except getopt.GetoptError:
	print('pick_samples_from_hybpiper.py -i <inputdir1/,inputdir2/,...> -o <outputdir/> -s <samplelist> -g <genelist>')
	sys.exit(2)
for opt, arg in opts:
	if opt == '-h':
		print('pick_samples_from_hybpiper.py -i <inputdir1/,inputdir2/,...> -o <outputdir/> -s <samplelist> -g <genelist>')
		sys.exit()
	elif opt in ("-i"):
		inputfiledir = arg.split(',')
	elif opt in ("-o"):
		outputfiledir = arg
	elif opt in ("-s"):
		samplelistfile = arg
	elif opt in ("-g"):
		genelistfile = arg

print('inputdirs: ')
print(inputfiledir)

# read file names from input directory

alignmentfilenames= []
for x in range(0,len(inputfiledir)):
	alignmentfilenames.append(glob.glob(inputfiledir[x]+'*.fasta'))

# read list of samples

with open(samplelistfile, 'r') as f:
	samplenames = f.read().splitlines()

print('selected samples:')
print(samplenames)

# read list of genes

with open(genelistfile, 'r') as f:
	genenames = f.read().splitlines()

print('selected genes:')
print(genenames)

# now that we have all the info we need:
# loop through input directories
#	loop through files found in input directories
#		if one contains a gene name, open it
#		loop through sequences in that file
#			if sequence contains one of the sample names, pick it

newalignments = [[] for n in range(len(genenames))]

# loop through directories
for dirx in range(0,len(alignmentfilenames)):
	# loop through files in directory
	for filex in range(0,len(alignmentfilenames[dirx])):
		# is this one of the gene files we want?
		if any(gn in alignmentfilenames[dirx][filex] for gn in genenames):
			# it is, but which gene? wonder if there is a more elegant way of doing this...
			whichgene = 0
			for genex in range(0,len(genenames)):
				if genenames[genex] in alignmentfilenames[dirx][filex]:
					whichgene = genex
			#now read sequence file and loop through sequences
			handle = open(alignmentfilenames[dirx][filex], "r")
			records = list(SeqIO.parse(handle,"fasta"))
			handle.close()
			for seqx in range(0,len(records)):
				if any(sn in records[seqx].name for sn in samplenames):
					newalignments[whichgene].append(records[seqx])

#write output alignments
for genex in range(0,len(genenames)):
	if len(newalignments[genex])>0:
		#print(newalignments[genex])
		SeqIO.write(newalignments[genex], outputfiledir+genenames[genex]+'.paralogs.fasta', "fasta")
