import sys, getopt
import glob
from Bio import AlignIO
from Bio import SeqIO

#list of variables to be populated from parameters
inputfiledir = '' #-i
outputfiledir = '' #-o

#parse parameters

try:
	opts, args = getopt.getopt(sys.argv[1:],"hi:o:")
except getopt.GetoptError:
	print('Script for reorganising PAFTOL sequences downloaded from the Kew Tree Explorer')
	print('from one file per sample into one file per gene, so that they have the same format')
	print('as HybPiper outputs. The input are files produced by process_paftol_outgroups_folder.py.')
	print('The outputs are for the sample picker script.')
	print('Use:')
	print('reorganise_paftol_seqs.py -i <inputdir1/,inputdir2/,...> -o <outputdir/>')
	sys.exit(2)
for opt, arg in opts:
	if opt == '-h':
		print('Script for reorganising PAFTOL sequences downloaded from the Kew Tree Explorer')
		print('from one file per sample into one file per gene, so that they have the same format')
		print('as HybPiper outputs. The input are files produced by process_paftol_outgroups_folder.py.')
		print('The outputs are for the sample picker script.')
		print('Use:')
		print('reorganise_paftol_seqs.py -i <inputdir1/,inputdir2/,...> -o <outputdir/>')
		sys.exit()
	elif opt in ("-i"):
		inputfiledir = arg.split(',')
	elif opt in ("-o"):
		outputfiledir = arg

print('inputdirs: ')
print(inputfiledir)

# read file names from input directory

alignmentfilenames= []
for x in range(0,len(inputfiledir)):
	alignmentfilenames.append(glob.glob(inputfiledir[x]+'*.fasta'))
print(alignmentfilenames)

# load all the alignments. I know that takes a lot of memory, but it seems best to get the gene list done first

alignments = []
genenames = set([])
samplenames = set([])
for dirx in alignmentfilenames:
	for filex in dirx:
		handle = open(filex, "r")
		records = list(SeqIO.parse(handle,"fasta"))
		handle.close()
		alignments.append(records)
		for record in records:
			print(record.name)
			genenames = genenames.union(set([record.name.split('-')[1]]))
			samplenames = samplenames.union(set([record.name.split('-')[0]]))
genenames = list(genenames)
samplenames = list(samplenames)

newalignments = [[] for n in range(len(genenames))]

for alignment in alignments:
	for record in alignment:
		whichgene = record.name.split('-')[1]
		samplename = record.name.split('-')[0]
		#while samplename[len(samplename)-1] == '.':
		#	samplename = samplename[:(len(samplename)-1)]
		samplename = samplename.replace('.','')
		whichgeneno = genenames.index(whichgene)
		record.name = ''
		record.id = samplename
		record.description = ''
		newalignments[whichgeneno].append(record)

#write output alignments
for genex in range(0,len(genenames)):
	if len(newalignments[genex])>0:
		#print(newalignments[genex])
		SeqIO.write(newalignments[genex], outputfiledir+genenames[genex]+'.fasta', "fasta")

samplenamefile = open(outputfiledir+'samplelist.txt', "w")
for samplename in samplenames:
	samplenamefile.write(samplename.replace('.','')+'\n')
samplenamefile.close()
