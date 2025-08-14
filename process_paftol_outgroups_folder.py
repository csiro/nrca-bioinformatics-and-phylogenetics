import string
from Bio import SeqIO
import sys
import glob


print('Script for reformatting PAFTOL sequences downloaded from Kew Tree Viewer.')
print('Input are folders containing fasta sequences for each sample with the name format')
print('INSDC.ERR0000000.Genus_species.a353.fasta')
print('Output are fasta files for each sample with more useful sequence names that can be.')
print('as input for the reorganise_paftol_seqs.py script.')
print('Command line parameter is input directory with .fasta files.')

pathtogenes = sys.argv[1]
if pathtogenes[len(pathtogenes)-1] != "/":
	pathtogenes = pathtogenes + "/"

targetsfilename = glob.glob(pathtogenes+'*.fasta')

print(targetsfilename)
samplenames = []

for i in range(0,len(targetsfilename)):
	handle = open(targetsfilename[i], "r")				#read input fasta file
	targets = list(SeqIO.parse(handle,"fasta"))
	handle.close()

	#genusnames=['']*len(records)					#create empty list for genus names of each sequence
	genenames=[]

	for x in range(0,len(targets)):						#make list of gene names by splitting sequence names in target file
		nameparts = targets[x].description.split(" ")
		taxon = nameparts[2][8:]
		if taxon[len(taxon)-1] == '.':
			taxon = taxon[:len(taxon)-1]
		targets[x].description = taxon+'_'+nameparts[4][12:]+'-'+nameparts[0]
		targets[x].id = ''
		targets[x].name = ''
		if (x==0):
			thisname = taxon+'_'+nameparts[4][12:]
			samplenames.append(thisname)

	filenameonly = targetsfilename[i].split('/')
	handle = open(targetsfilename[i][:(len(targetsfilename[i])-len(filenameonly[len(filenameonly)-1]))] + "edited_"+ filenameonly[len(filenameonly)-1],"w")

	#for x in range(0,len(targets)):          		        #loop through seqs
	SeqIO.write(targets, handle, "fasta")			#add target sequence to alignment
	handle.close()

handle2 = open(pathtogenes+"samples.txt","w")
for i in range(0,len(samplenames)):
	handle2.write(samplenames[i]+'\n')
handle2.close()
