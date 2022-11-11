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

# get all fasta files in folder
fastafilename = glob.glob(pathtogenes+"*"+fastaext)

# run iqtree for all alignments
for i in range(0, len(fastafilename)):
	#subprocess.run(["mafft","--auto","--adjustdirectionaccurately",fastafilename[i],">",fastafilename[i][:(len(fastafilename[i])-len(fastaext))]+".aligned"+fastaext])
	mycmd = "mafft --auto --adjustdirection "+fastafilename[i]+" > "+fastafilename[i][:(len(fastafilename[i])-len(fastaext))]+".aligned"+fastaext
	subprocess.run(mycmd,shell=True)
