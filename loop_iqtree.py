import subprocess
import glob
import sys

# look at command line parameters
if len(sys.argv)>1:
	pathtogenes = sys.argv[1]		#first command line parameter may be path to folder with alignments, but is optional
else:
	pathtogenes = "./"				#if not provided, it is assumed that the folder to be worked in is .
if pathtogenes[len(pathtogenes)-1] != "/":
	pathtogenes = pathtogenes + "/"

if len(sys.argv)>2:
	fastaext = sys.argv[2]			# second command line parameter is extension of fasta files, but is optional
else:
	fastaext = ".fasta"				# if not provided, extension is assumed to be .fasta

if len(sys.argv)>3:
	if sys.argv[3] == 'y':
		parts = True			# flag whether to look for partition file with same name as input file except partition instead of fasta at end
	else:
		parts = False
else:
	parts = False				# if not provided, assume there is no such file

if len(sys.argv)>4:
	otherparams = sys.argv[4]		# any other parameters for IQ-TREE
else:
	otherparams = ""

# get all fasta files in folder
fastafilename = glob.glob(pathtogenes+"*"+fastaext)

# run iqtree for all alignments
for i in range(0, len(fastafilename)):
	if parts:
		mycmd = 'iqtree2 -s '+ fastafilename[i] + ' -spp ' + fastafilename[i][:(len(fastafilename[i])-6)] + '.partition '+ otherparams
	else:
		mycmd = 'iqtree2 -s '+ fastafilename[i] + ' '+ otherparams
	subprocess.call(mycmd,shell=True)
