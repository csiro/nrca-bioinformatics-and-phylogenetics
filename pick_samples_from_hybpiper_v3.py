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
	print('pick_samples_from_hybpiper.py -i <inputdirlist> -o <outputdir/> -s <sample whitelist> -g <genelist>')
	print('This version expects a list of input directories in a text file. If it finds the same sample in several directories, it will ignore all instances except the first.')
	sys.exit(2)
for opt, arg in opts:
	if opt == '-h':
		print('pick_samples_from_hybpiper.py -i <inputdirlist> -o <outputdir/> -s <sample whitelist> -g <genelist>')
		print('This version expects a list of input directories in a text file. If it finds the same sample in several directories, it will ignore all instances except the first.')
		sys.exit()
	elif opt in ("-i"):
		inputdirfile = arg
	elif opt in ("-o"):
		outputfiledir = arg
	elif opt in ("-s"):
		samplelistfile = arg
	elif opt in ("-g"):
		genelistfile = arg

# read list of samples

with open(inputdirfile, 'r') as f:
	inputfiledir = f.read().splitlines()

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

# prepare list to check if sample has actually been found

samplecount = [0]*len(samplenames)
samplesfoundglobal = set([])

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
	samplesfoundinthisdir = set([])
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
				currentsamplename = records[seqx].name.split('.')[0].split(' ')[0]
				#if any(sn in records[seqx].name for sn not in samplenames):
				if (currentsamplename in samplenames):
					samplesfoundinthisdir.add(currentsamplename)
					if (samplecount[samplenames.index(currentsamplename)] == 0):
						newalignments[whichgene].append(records[seqx])
	#samplesfoundinthisdir = set(samplesfoundinthisdir)

	#print("samples in "+inputfiledir[dirx])
	#print(samplesfoundinthisdir)

	samplesfoundglobal.update(samplesfoundinthisdir)

	#print("global samples found")
	#print(samplesfoundglobal)
	
	for samplesx in range(0,len(samplenames)):
		if (samplenames[samplesx] in samplesfoundinthisdir):
			samplecount[samplesx] = samplecount[samplesx] + 1

#write output alignments
for genex in range(0,len(genenames)):
	if len(newalignments[genex])>0:
		#print(newalignments[genex])
		SeqIO.write(newalignments[genex], outputfiledir+genenames[genex]+'.paralogs.fasta', "fasta")

#alert user to samples found or samples found repeatedly

print('\nThe following samples were not found:')
for samplex in range(0, len(samplenames)):
	if samplecount[samplex] == 0:
		print(samplenames[samplex])
print('\nThe following samples were found in several folders:')
for samplex in range(0, len(samplenames)):
	if samplecount[samplex] > 1:
		print(samplenames[samplex])

#print('')
#for samplex in range(0, len(samplenames)):
#	print(samplenames[samplex]+'   '+str(samplecount[samplex]))
