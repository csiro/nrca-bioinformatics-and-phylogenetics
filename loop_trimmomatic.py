import subprocess
import glob
import sys

# look at command line parameters
pathtogenes = sys.argv[1]		#first command line parameter is path to folder with read files
if pathtogenes[len(pathtogenes)-1] != "/":
	pathtogenes = pathtogenes + "/"

if len(sys.argv)>2:
	otherparams = sys.argv[2]		# any other parameters for Trimmomatic
else:
	otherparams = "CROP:147 HEADCROP:3 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50"

# get all fasta files in folder
r1files = glob.glob(pathtogenes+"*_R1.fastq")
r2files = glob.glob(pathtogenes+"*_R2.fastq")
print(r1files)
r1files.sort()
r2files.sort()
print(r1files)

# run iqtree for all alignments
for i in range(0, len(r1files)):
	r1outfilename = r1files[i][:(len(r1files[i])-6)]
	r2outfilename = r2files[i][:(len(r2files[i])-6)]
	mycmd = 'trimmomatic PE '+r1files[i]+' '+r2files[i]+' '+r1outfilename+'_FP.fastq '+r2outfilename+'_FU.fastq '+r1outfilename+'_RP.fastq '+r2outfilename+'_RU.fastq '+otherparams
	subprocess.call(mycmd,shell=True)
