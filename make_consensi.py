from Bio import AlignIO
from Bio import SeqIO
import string, sys, getopt, glob, os
from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC
import Bio.Align
from Bio.SeqRecord import SeqRecord
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment

# defaults
pathtogenes = './' #-i
fastaextension = '.fasta' #-e
outputfile = "./consensi.fasta" #-o
seqnametxt = ''

#parse parameters

try:
	opts, args = getopt.getopt(sys.argv[1:],"hi:o:e:s:")
except getopt.GetoptError:
	print('make_consensi.py -i <inputdir/> -o <outputfile> -e <inputfileextension> -s <species>')
	print('More details: make_consensi.py -h')
	sys.exit(2)
for opt, arg in opts:
	if opt == '-h':
		print('Script for creating consensus sequences for a species or an entire alignment, across a folder containing multiple fasta files.')
		print('\nUse: make_consensi.py -i <inputdir/> -o <outputfile> -e <inputfileextension> -s <species>')
		print('\nParameters')
		print('e: Extension of input (fasta) files to be read, default is .fasta')
		print('h: Help, i.e. get this information')
		print('i: Path to directory containing input (fasta) files, default is current working directory')
		print('o: Output file; default is consensi.fasta in current working directory')
		print('s: Text that is part of sequence names to decide which sequences the consensus sequences will be created from. If left empty, the consensus will be created from all sequences in the alignment.')
		sys.exit()
	elif opt in ("-i"):
		pathtogenes = arg
		if pathtogenes[len(pathtogenes)-1] != "/":
			pathtogenes = pathtogenes + "/"
	elif opt in ("-o"):
		outputfile = arg
	elif opt in ("-e"):
		fastaextension = arg
	elif opt in ("-s"):
		seqnametxt = arg

# glob all the input files assuming they end in fasta
inputfilenames = glob.glob(pathtogenes + "*"+fastaextension)

all_consensi = []

for i in range(0,len(inputfilenames)):
	# read alignment
	alignmentfilename = inputfilenames[i]
	alignment = AlignIO.read(alignmentfilename, "fasta")
	#filter out anything that doesn't contain the specified sequence name part
	if seqnametxt != '':
		culled_align = Bio.Align.MultipleSeqAlignment([])
		for j in range(0,len(alignment)):
			if seqnametxt in alignment[j].id:
				culled_align.append(alignment[j])
	else:
		culled_align = alignment
	if len(culled_align) > 0:
		summary_aln = AlignInfo.SummaryInfo(culled_align)
		myconsseq = summary_aln.gap_consensus(threshold=0.4,ambiguous='n')
		mycons = SeqRecord(myconsseq, id=os.path.basename(inputfilenames[i])[0:4], name=os.path.basename(inputfilenames[i])[0:4],description="")
		all_consensi.append(mycons)
Bio.SeqIO.write(all_consensi,outputfile,format='fasta')


