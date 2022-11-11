from Bio import AlignIO
from Bio import SeqIO
import sys, getopt, subprocess
import string
import glob
from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC
import Bio.Align
from Bio.SeqRecord import SeqRecord

# function to translate DNA into AA, turning partial codons into "?"
def mytranslate(aseqrecord):
	codons = ["TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG", "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG", "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG", "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG", "---"]
	AAs = ["F", "F", "L", "L", "S", "S", "S", "S", "Y", "Y", "X", "X", "C", "C", "X", "W", "L", "L", "L", "L", "P", "P", "P", "P", "H", "H", "Q", "Q", "R", "R", "R", "R", "I", "I", "I", "M", "T", "T", "T", "T", "N", "N", "K", "K", "S", "S", "R", "R", "V", "V", "V", "V", "A", "A", "A", "A", "D", "D", "E", "E", "G", "G", "G", "G", "-"]
	recaalength = int(len(aseqrecord.seq) / 3)
	atranslation = ""
	for i in range(0,recaalength):
		if (i*3)+2 > len(aseqrecord.seq):
			atranslation = atranslation + "?"
		else:
			thiscodon = aseqrecord.seq[(i*3):((i+1)*3)]
			thiscodon = thiscodon.upper()
			if thiscodon in codons:
				atranslation = atranslation + AAs[codons.index(thiscodon)]
			else:
				atranslation = atranslation + "?"
	return SeqRecord(Seq(atranslation),id=aseqrecord.id,description='')

# function to count number of variable characters (columns) in an alignment
def count_variable_chars(seqalignment):
	diffcount = 0
	for x in range(0,seqalignment.get_alignment_length()):
		columndata = set()
		for y in range(0,len(seqalignment)):
			columndata.add(seqalignment[y].seq[x])
		# remove ambiguous or missing data so they don't get counted as variation
		columndata.discard('-')
		columndata.discard('X')
		columndata.discard('?')
		if len(columndata)>1:
			diffcount = diffcount + 1
	return diffcount

#list of variables to be populated from parameters, with defaults, excepting the targets file, which doesn't have one
pathtogenes = './' #-i
outputfiledir = './' #-o
fastaextension = '.fasta' # -e
mingaplength = 15 # -g
matchvalue = 5 # -m

#parse parameters

try:
	opts, args = getopt.getopt(sys.argv[1:],"hi:o:e:r:m:g:")
except getopt.GetoptError:
	print('check_odd_indels_and_frames2.py -i <inputdir/> -o <outputdir/> -e <inputfileextension> -r <referencefile> -g <mingaplength> -m <matchvalue>')
	print('More details: check_odd_indels_and_frames2.py -h')
	sys.exit(2)
for opt, arg in opts:
	if opt == '-h':
		print('Script for ensuring that multiple gene alignments are cleaned of dirty assembly edges and in frame before concatenation. Requires input folder with fasta gene alignment files and a fasta file with a guaranteed to be in frame, trusted reference sequence for each gene. The names of the latter should be formatted as targetname-geneidentifier, and the gene identifier must be part of the relevant gene alignment file in the input directory. Also requires mafft to be installed and in the path, to allow alignment.')
		print('\nUse: check_odd_indels_and_frames2.py -i <inputdir/> -o <outputdir/> -e <inputfileextension> -r <referencefile> -g <mingaplength> -m <matchvalue>')
		print('\nParameters')
		print('e: Extension of input (fasta) files to be read, default is .fasta')
		print('g: Minimum length of gap (indel) for its edges to be examined for match to reference, default is 15.')
		print('h: Help, i.e. get this information')
		print('i: Path to directory containing input (fasta) files, default is current working directory')
		print('m: Number of bases that need to match reference around a gap or at ends for cutting to stop, default is 5')
		print('o: Path to output directory, default is current working directory')
		print('r: Name of file containing reference sequences; this argument is required, no default')
		sys.exit()
	elif opt in ("-i"):
		pathtogenes = arg
		if pathtogenes[len(pathtogenes)-1] != "/":
			pathtogenes = pathtogenes + "/"
	elif opt in ("-o"):
		outputfiledir = arg
		if outputfiledir[len(outputfiledir)-1] != "/":
			outputfiledir = outputfiledir + "/"
	elif opt in ("-e"):
		fastaextension = arg
	elif opt in ("-r"):
		targetsfilename = arg

# glob all the input files assuming they end in fasta
inputfilenames = glob.glob(pathtogenes + "*"+fastaextension)

# add references to all the genes

handle = open(targetsfilename, "r")				#read fasta file with references
targets = list(SeqIO.parse(handle,"fasta"))
handle.close()

#genusnames=['']*len(records)					#create empty list for sample names of each sequence
genenames=[]

for x in range(0,len(targets)):						#make list of gene names by splitting sequence names in target file
	nameparts = targets[x].description.split("-")
	genenames.append(nameparts[1])

#for x in range(0,len(genenames)):          		        #loop through genes
#	if os.path.isfile(pathtogenes+genenames[x]+fasta_name_end):
#		handle = open(pathtogenes+genenames[x]+fasta_name_end,"a")		#open gene alignment for appending; text at end either ".paralogs.fasta" or ".FNA"
#		# added the name change to remove gene name - check if that causes problems later, but I guess it shouldn't
#		targets[x].id = targets[x].id.split('-')[0]
#		targets[x].name = ''
#		targets[x].description = ''
#		SeqIO.write(targets[x], handle, "fasta")			#add target sequence to alignment
#		handle.close()


reportfile = open(outputfiledir+'check_odd_indels_and_frames_report.txt','w')

for i in range(0,len(inputfilenames)):
	# read alignment
	alignmentfilename = inputfilenames[i]
	reportfile.write('Reading '+alignmentfilename+'\n')
	alignment = AlignIO.read(alignmentfilename, "fasta")

	# find out which reference to use for this gene
	target = -99
	for x in range(0,len(genenames)):
		if genenames[x] in alignmentfilename:
			target = x
	if target != -99:
		# add reference sequence
		# SHITE, PROBLEM IS IF THE TARGET IS LONGER THAN THE REFERENCE. MAYBE CHANGE ALL TO SEQ INSTEAD OF ALIGN...
		if len(targets[target].seq) < len(alignment[0].seq):
			targets[target].seq = targets[target].seq + ('-'*(len(alignment[0].seq)-len(targets[target].seq)))
		if len(targets[target].seq) > len(alignment[0].seq):
			for y in range(0,len(alignment)):
				alignment[y].seq = alignment[y].seq + ('-'*(len(targets[target].seq)-len(alignment[0].seq)))
		alignment.append(targets[target])
		target = len(alignment)-1
		
		# realign: save alignment, call MAFFT, load result
		with open(outputfiledir+'temporary.fasta', "w") as output_handle:
			SeqIO.write(alignment, output_handle, "fasta")
		mycmd = "mafft --auto --adjustdirection --inputorder "+outputfiledir+"temporary.fasta > "+outputfiledir+"temporary.aligned.fasta"
		subprocess.run(mycmd,shell=True)
		alignment = AlignIO.read(outputfiledir+'temporary.aligned.fasta', "fasta")

		# hm. turns out mafft adds "R_" to sequences it had to reverse-complement
		for x in range(0,len(alignment)):
			if alignment[x].id[:3] == "_R_":
				alignment[x].name = alignment[x].id[3:]
				alignment[x].name=''
				alignment[x].description=''
		
		# trim alignment to target, from left and right sides
		x = 0
		while alignment[target].seq[x] =="-":
			x = x + 1
		if x > 0:
			alignment = alignment[:,x:]
		x = alignment.get_alignment_length()-1
		while alignment[target].seq[x] =="-":
			x = x - 1
		if x < alignment.get_alignment_length()-1:
			alignment = alignment[:,:x]
		# walk through target sequence and see where it has gaps that cannot be divided by three
		# then slice out those gaps, so that we can partition the analysis by codon position
		currentlength = alignment.get_alignment_length() 
		y = 0
		gapcounter = 0
		while y < currentlength:
			if alignment[target].seq[y] == "-":
				gapcounter = gapcounter + 1
				y = y + 1
			else:
				if gapcounter % 3 != 0:
					# slice out offending gap, recalculate alignment length, adjust y
					for x in range(0,len(alignment)):
						alignment[x].seq =  alignment[x].seq[:(y-gapcounter)] + alignment[x].seq[y:]
					y = y - gapcounter
					currentlength = alignment.get_alignment_length()
				else:
					y = y + 1
				gapcounter = 0
		
		# NEED TO ADD THE ALGORITHM THAT REMOVES BAD ENDS AROUND A GAP !!!!!!!!!!!!!!!!!1
		# begin cleaning bad ends
		
		didicut = False
		for c in range(0,len(alignment)):
			if c != target:
				# walk through alignment from the left
				# ensure that cuts will be made directly from beginning of alignment and after longer gaps,
				# but not after small gaps
				x = 0
				longgap = True
				while x < (len(alignment[c].seq)-1):
					gaplengthcount = 0
					while (alignment[c].seq[x] =="-") & (x < (len(alignment[c].seq)-1)):
						x = x + 1
						gaplengthcount = gaplengthcount + 1
						#print("counter: "+str(x))
					if gaplengthcount > mingaplength:		# change this to adjust minimum gap length for cutting at the edges
						longgap = True
					if longgap & (x < (len(alignment[c].seq)-1-matchvalue)):
						# keep turning bases into gaps until at least "matchvalue" consecutive ones were identical to target
						shouldicut = False
						startingmatch = x
						x = x + matchvalue
						while (x < len(alignment[c].seq)-1) & (alignment[c].seq[(x-matchvalue):x].upper() != alignment[target].seq[(x-matchvalue):x].upper()):
							x = x + 1
							shouldicut = True
						if shouldicut:
							alignment[c].seq = alignment[c].seq[:startingmatch]+'-'*(x-startingmatch-matchvalue)+alignment[c].seq[(x-matchvalue):]
							dididcut = True
							#for y in range(startingmatch,x):
							#	alignment[i].seq[:y] = "-"
					# move to end of remaining sequence section
					x = x + 1
					if x < (len(alignment[c].seq)-1):
						while (alignment[c].seq[x] != "-") & (x < (len(alignment[c].seq)-2)):
							x = x + 1
					longgap = False
				# now the same from the right end
				x = alignment.get_alignment_length()-1
				longgap = True
				while x > 0:
					gaplengthcount = 0
					while (alignment[c].seq[x] =="-") & (x > 1):
						x = x - 1
						#print("counter: "+str(x))
						gaplengthcount = gaplengthcount + 1
					if gaplengthcount > mingaplength:		# change this to adjust minimum gap length for cutting at the edges
						longgap = True
					if longgap & (x > matchvalue):
						# keep turning bases into gaps until at least "matchvalue" consecutive ones were identical to target
						shouldicut = False
						startingmatch = x
						x = x - matchvalue
						while alignment[c].seq[x:(x+matchvalue)].upper() != alignment[target].seq[x:(x+matchvalue)].upper():
							x = x - 1
							shouldicut = True
						if shouldicut:
							alignment[c].seq = alignment[c].seq[:(x+matchvalue)]+'-'*(startingmatch-x-matchvalue+1)+alignment[c].seq[(startingmatch+1):]
							dididcut = True
							#for y in range(x, (startingmatch+1)):
							#	alignment[c].seq[y] = "-"
					# move to end of remaining sequence section
					if x > 0:
						x = x - 1
						while (alignment[c].seq[x] != "-") & (x > 1):
							x = x - 1
					longgap = False
		if didicut:
			reportfile.write("Removed bad ends.\n")

		# end cleaning up bad ends

	else:
		reportfile.write('WARNING: no reference for this gene file\n')

	# check all three frames for their number of variable characters
	mostpars = 0
	parsscore = 999999
	for framecount in range(0,3):
		aa_align = Bio.Align.MultipleSeqAlignment([])
		for i in range(0,len(alignment)):
			thisseq = SeqRecord(alignment[i].seq[framecount:], alignment[i].id)
			# attach translation to growing translated alignment, then check for differences
			aa_align.append(mytranslate(thisseq))
		thisparsscore = count_variable_chars(aa_align)
		if thisparsscore < parsscore:
			parsscore = thisparsscore
			mostpars = framecount
	# unfortunately I have no idea how to add a column of Ns, so the easiest way to shift the frame
	# is to delete columns at the beginning
	alignment = alignment[:,mostpars:]
	# now delete columns at the end to allow simpler codon partitioning in concatenated DNA sequence alignment
	alignment = alignment[:,:(int(alignment.get_alignment_length()/3)*3)]
	# get again the translation with the lowest number of variable characters
	aa_align = Bio.Align.MultipleSeqAlignment([])
	for i in range(0,len(alignment)):
		thisseq = SeqRecord(alignment[i].seq, id=alignment[i].id, description=alignment[i].description)
		aa_align.append(mytranslate(thisseq))
	# remove last sequence, because that is the reference
	if target != -99:
		alignment = alignment[:(len(alignment)-1),:]
		aa_align = aa_align[:(len(aa_align)-1),:]
	
	# remove gap-only columns, at least for entire codons (I think there shouldn't be any 
	currentpos = 0
	while currentpos < aa_align.get_alignment_length():
		thiscol = list(set(aa_align[:,currentpos]))
		if ((len(thiscol)==1) and (thiscol[0]=='-')):
			alignment = alignment[:,:(3*currentpos)] + alignment[:,((3*currentpos)+3):]
			aa_align = aa_align[:,:currentpos] + aa_align[:,(currentpos+1):]
		else:
			currentpos = currentpos + 1

	# write alignment
	nameholder = alignmentfilename.split('/')
	alignmentfilenameonly = nameholder[len(nameholder)-1]
	AlignIO.write(alignment, outputfiledir+alignmentfilenameonly[:len(alignmentfilenameonly)-len(fastaextension)]+".frame"+fastaextension, "fasta")
	AlignIO.write(aa_align, outputfiledir+alignmentfilenameonly[:len(alignmentfilenameonly)-len(fastaextension)]+".AA"+fastaextension, "fasta")
reportfile.close()
