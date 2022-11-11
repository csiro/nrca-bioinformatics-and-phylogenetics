import string
import glob
from Bio import SeqIO
from Bio import Seq
import sys, getopt

#defaults for parameters

pathtogenes = "./"   # -i
fastaextension = ".fasta"   # -e
outputfiledir = "./"  # -o
outfilename = "concat"  # -n
nameseparator = "_"  # -u
includesamples =[]  # (-s provides filename for list)
mingenes = 1  # -g
mintaxa = 5  # -t
identifyingnameparts = [0,1]  # -c

#parse parameters

try:
	opts, args = getopt.getopt(sys.argv[1:],"hc:e:i:g:n:o:s:t:u:")
except getopt.GetoptError:
	print('Use: concat7.py -i <inputdirectory> -o <outputdirectory> -e <inpufileextension> -n <outputfilename> -u <nameseparator> -c <identifyingnameparts> -s <samplelist> -t <mintaxa> -g <mingenes>')
	print('More details: concat7.py -h')
	sys.exit(2)
for opt, arg in opts:
	if opt == '-h':
		print('Script for concatenation of sequence files for phylognenetic analysis. Input is a folder with fasta files in which the samples must have the same name across alignments, and optionally a text file containing the list of samples to include. Outputs are concatenated alignments in PHYLIP, Nexus and TNT formats, partition files, gene alignments for only the selected files plus codon partition files, sample assignment information for ASTRAL, and a table showing the numbers of genes for each sample and vice versa.')
		print('\nUse: concat7.py -i <inputdirectory> -o <outputdirectory> -e <inpufileextension> -n <outputfilename> -u <nameseparator> -c <identifyingnameparts> -s <samplelist> -t <mintaxa> -g <mingenes>')
		print('\nParameters')
		print('c: Position of name parts that identify species for sample-to-species lumping for ASTRAL, starting count with zero. E.g. for Genus_species_sampleID, it would be 0,1 (note absence of whitespace). Always requires two numbers, use same number if only one is required. Default is 0,1')
		print('e: Extension of input (fasta) files to be read, default is .fasta')
		print('i: Path to directory containing input (fasta) files, default is current working directory')
		print('g: Minimum number of genes for sample to be included, default is 1, but recommended to be more restrictive')
		print('h: Help, i.e. get this information')
		print('n: Beginning of output file names, default is concat')
		print('o: Path to output directory, default is current working directory')
		print('s: Name of (including path to) text file with names of samples to include, default is none (to pick all samples)')
		print('t: Minimum number of taxa for gene to be included, default is 5')
		print('u: Separator of name parts in sample names, default is underscore as in Genus_species_sampleID\n')
		sys.exit()
	elif opt == "-i":
		pathtogenes = arg
		if pathtogenes[len(pathtogenes)-1] != "/":
			pathtogenes = pathtogenes + "/"
	elif opt in ("-e"):
		fastaextension = arg
	elif opt in ("-o"):
		outputfiledir = arg
		if outputfiledir[len(outputfiledir)-1] != "/":
			outputfiledir = outputfiledir + "/"
	elif opt in ("-n"):
		outfilename = arg
	elif opt in ("-u"):
		nameseparator = arg
	elif opt in ("-s"):
		includefilename = arg
		with open(includefilename, 'r') as f:
			includesamples = f.read().splitlines()
	elif opt in ("-g"):
		mingenes = int(arg)
	elif opt in ("-t"):
		mintaxa = int(arg)
	elif opt in ("-c"):
		identifyingnameparts = arg.split(',')[:2]
		identifyingnameparts = [int(a) for a in identifyingnameparts]

print('pathtogenes '+pathtogenes)
print('extension '+fastaextension)

#query names of files to be concatenated

nextname ="dummy"
inputfilenames = []

inputfilenames = glob.glob(pathtogenes+"*"+fastaextension)

print(inputfilenames)	
print(len(inputfilenames))

#read files

datasets = []

for x in range(0,len(inputfilenames)):
	handle = open(inputfilenames[x], "r")				#read input fasta file
	records = list(SeqIO.parse(handle,"fasta"))
	datasets += [records]
	handle.close()

# ensure that leading _R_ is taken off, just in case
for x in range(0,len(datasets)):
	for y in range(0,len(datasets[x])):
		if datasets[x][y].id[:3] == "_R_":
			datasets[x][y].id = datasets[x][y].id[3:]
		if datasets[x][y].id[:3] == "_R_":
			datasets[x][y].id = datasets[x][y].id[3:]

# get all sample names across all genes, then reduce to uniques for concatenation
allsamplenames = []
samplesbygene = []
for x in range(0,len(inputfilenames)):
	nameholder = []
	for y in range (0, len(datasets[x])):
		nameholder += [datasets[x][y].id]
		allsamplenames += [datasets[x][y].id]
	samplesbygene.append(set(nameholder))

allsamplenames = list(set(allsamplenames))

# now drop excluded terminals from list

if len(includesamples) > 0:
	excludesamples = []
	for x in range(0,len(allsamplenames)):
		if allsamplenames[x] not in includesamples:
			excludesamples.append(allsamplenames[x])
	if len(excludesamples) > 0:
		for x in range(0,len(excludesamples)):
			allsamplenames.remove(excludesamples[x])

# iteratively reduce data matrix to ensure minimum coverage across genes and terminals

print('Mintaxa: '+str(mintaxa))
print('Mingenes: '+str(mingenes))

genesremoved = []
taxaremoved = []
happynow = False
while not(happynow):
	happynow = True
	taxcoverage = [0]*len(allsamplenames)
	genecoverage = [0]*len(datasets)
	for x in range(0, len(datasets)):
		for y in range(0, len(allsamplenames)):
			if allsamplenames[y] in samplesbygene[x]:
				taxcoverage[y] += 1
				genecoverage[x] += 1
	#print('Genecoverage:\n')
	#print(genecoverage)
	#print('\nTaxon coverage:\n')
	#print(taxcoverage)
	x = 0
	while x < (len(datasets)):
		if (genecoverage[x] < mintaxa):
			happynow = False
			datasets.pop(x)
			genecoverage.pop(x)
			samplesbygene.pop(x)
			genesremoved.append(inputfilenames.pop(x))
		else:
			x += 1
	x = 0
	while x < (len(allsamplenames)):
		if (taxcoverage[x] < mingenes):
			happynow = False
			taxaremoved.append(allsamplenames.pop(x))
			taxcoverage.pop(x)
		else:
			x += 1

#figure out how long the sequences in each dataset are - assuming they are all equally long in each of them!

nchar = 0
seqlengths = [0]*len(inputfilenames)
for d in range(0,len(inputfilenames)):
	seqlengths[d]=len(datasets[d][0].seq)
	nchar = nchar + seqlengths[d]
	for x in range(0, len(datasets[d])):
		if (len(datasets[d][x].seq) != seqlengths[d]):
			print("Error: "+datasets[d][x].name+" in "+inputfilenames[d] + " has deviating sequence length. Added question marks to end of sequence.\n")
			corrseq = ''.join(list(datasets[d][x].seq))
			corrseq = corrseq + ''.join('?'*int(seqlengths[d]-len(datasets[d][x].seq)))
			datasets[d][x].seq = Seq.Seq(corrseq)

#finaldata = [["?"]*len(inputfilenames)]*len(allsamplenames)

#compose data matrix

finaldata = []

for x in range(0,len(allsamplenames)):     #loop through taxa
	currentdataline = []
	for d in range(0, len(inputfilenames)):    #loop through datasets
		foundflag=0
		for y in range(0, len(datasets[d])):      #find taxon in dataset
			if (datasets[d][y].name == allsamplenames[x]):
				#finaldata[x][d] = str(datasets[d][y].seq)
				currentdataline.append(str(datasets[d][y].seq))
				foundflag=1
		if foundflag == 0:
			currentdataline.append(''.join('-'*seqlengths[d]))
	finaldata.append(currentdataline)

#start writing output files

raxmloutfilename = outputfiledir + outfilename + ".phy"
raxmloutfile = open(raxmloutfilename, "w")
raxmloutfile.write(str(len(allsamplenames))+" "+str(nchar)+"\n")

partsoutfilename = outputfiledir + outfilename + "_raxml.parts"
partsoutfile = open(partsoutfilename, "w")

iqpartsoutfilename = outputfiledir + outfilename + "_iqtgenes.parts"
iqpartsoutfile = open(iqpartsoutfilename, "w")

nexusoutfilename = outputfiledir + outfilename + ".nex"
nexusoutfile = open(nexusoutfilename, "w")
nexusoutfile.write("#NEXUS\n")
nexusoutfile.write("BEGIN DATA;\n")
nexusoutfile.write("\tDIMENSIONS NTAX="+str(len(allsamplenames))+" NCHAR="+str(nchar)+";\n")
nexusoutfile.write("\tformat datatype=dna missing=? gap=-;\n")
nexusoutfile.write("MATRIX\n")

tntoutfilename = outputfiledir + outfilename + ".tnt"
tntoutfile = open(tntoutfilename, "w")
tntoutfile.write("xread\n")
tntoutfile.write(str(nchar)+" "+str(len(allsamplenames))+"\n")
tntoutfile.write("&[dna]\n")

fastaoutfilename = outputfiledir + outfilename + ".fasta"
fastaoutfile = open(fastaoutfilename, 'w')

#write data matrices

for x in range(0,len(allsamplenames)):
	raxmloutfile.write(allsamplenames[x]+"     ")
	nexusoutfile.write("\t"+allsamplenames[x]+"\t")
	tntoutfile.write(allsamplenames[x]+"\t")
	fastaoutfile.write('>'+allsamplenames[x]+'\n')
	raxmloutfile.write(''.join(finaldata[x]) + "\n")
	nexusoutfile.write(''.join(finaldata[x]) + "\n")
	tntoutfile.write(''.join(finaldata[x]) + "\n")
	fastaoutfile.write(''.join(finaldata[x]) +'\n')

#write partitions and gene alignments (the latter necessary because taxa and genes have been excluded)

codonpartsfile = open(outputfiledir + outfilename+'_iqtcodons.parts','w')    # write codon partition file for IQ-TREE
codonpartsfile.write('#nexus\nbegin sets;\n\tcharset part1 = 1-' +str(nchar)+ '\\3;\n\tcharset part2 = 2-' +str(nchar)+ '\\3;\n\tcharset part3 = 3-' +str(nchar)+ '\\3;\nend;\n')
codonpartsfile.close()

nexusoutfile.write(";\n")
nexusoutfile.write("END;\n")
nexusoutfile.write("begin sets;\n")

iqpartsoutfile.write('#nexus\nbegin sets;\n')

allsets = ''

currentsum = 0
for d in range(0, len(inputfilenames)):
	partsoutfile.write("DNA, gene"+str(d+1)+" = "+ str(currentsum+1) +"-"+ str(currentsum+seqlengths[d])+"\n")
	nameholder = inputfilenames[d][:(len(inputfilenames[d])-len(fastaextension))]    # remove extension
	nameholder = nameholder.split('/') # get rid of path
	nameholder = nameholder[len(nameholder)-1]
	namefornexus = nameholder.translate(str.maketrans('-. /\\','_____')) # get rid of name parts that might trouble phylogenetics software
	nexusoutfile.write("\tcharset "+namefornexus+" = "+ str(currentsum+1) +"-"+ str(currentsum+seqlengths[d])+";\n")
	iqpartsoutfile.write("\tcharset "+namefornexus+" = "+ str(currentsum+1) +"-"+ str(currentsum+seqlengths[d])+";\n")
	allsets = allsets + namefornexus +':'+str(currentsum+1) +"-"+ str(currentsum+seqlengths[d]) +','
	currentsum = currentsum + seqlengths[d]
	datasetholder = []
	for e in range(0,len(datasets[d])):
		if datasets[d][e].id in allsamplenames:
			datasetholder += [datasets[d][e]]
	geneoutfilename = inputfilenames[d][:(len(inputfilenames[d])-len(fastaextension))].split('/')
	geneoutfilename = outputfiledir+geneoutfilename[len(geneoutfilename)-1]
	SeqIO.write(datasetholder, geneoutfilename + '.leftover.fasta', "fasta")
	genepartoutfile = open(geneoutfilename + '.leftover.partition', 'w')
	genepartoutfile.write('#nexus\nbegin sets;\n\tcharset part1 = 1-' +str(len(datasetholder[0].seq))+ '\\3;\n\tcharset part2 = 2-' +str(len(datasetholder[0].seq))+ '\\3;\n\tcharset part3 = 3-' +str(len(datasetholder[0].seq))+ '\\3;\nend;\n')
	genepartoutfile.close()

nexusoutfile.write('\n\tcharpartition combined = '+allsets[:(len(allsets)-1)]+';\n')    # write combined gene partition for SVDq in PAUP
iqpartsoutfile.write('end;\n')

# now also species assignment of samples

astralfile = open(outputfiledir + outfilename+'_astralspecies','w')
allspeciesnames = []
for i in range(0,len(allsamplenames)):           # figure out unique species for assigning samples to species, assuming first two name parts are that information
	namesplit = allsamplenames[i].split(nameseparator)
	if len(namesplit) == 1:
		allspeciesnames += [namesplit[0]]
	else:
		allspeciesnames += [namesplit[identifyingnameparts[0]]+nameseparator+namesplit[identifyingnameparts[1]]]
allspeciesnames = list(set(allspeciesnames))
nexusoutfile.write('\ttaxpartition species = ')
for i in range(0,len(allspeciesnames)):
	if i>0:
		nexusoutfile.write(', ')
	nexusoutfile.write(allspeciesnames[i]+': ')
	astralfile.write(allspeciesnames[i]+':')
	commaflag = False
	for j in range(0,len(allsamplenames)):
		if allsamplenames[j][:len(allspeciesnames[i])] ==  allspeciesnames[i]:
			if commaflag:
				astralfile.write(',')
			else:
				commaflag = True
			nexusoutfile.write(allsamplenames[j]+' ')
			astralfile.write(allsamplenames[j])
	astralfile.write('\n')
nexusoutfile.write(';\n')
astralfile.close()

#finish writing outfiles

nexusoutfile.write("END;\n")
tntoutfile.write(";\n")

#close all outfiles

nexusoutfile.close()
raxmloutfile.close()
tntoutfile.close()
partsoutfile.close()
iqpartsoutfile.close()
fastaoutfile.close()

# write coverage stats

coveragefilename = outputfiledir + outfilename + "_coverage.tsv"
coveragefile = open(coveragefilename, "w")
coveragefile.write('Genes per terminal:\n')
for x in range(0,len(allsamplenames)):
	coveragefile.write(allsamplenames[x] + '\t' + str(taxcoverage[x]) + '\n')
coveragefile.write('\nTerminals removed because less than '+str(mingenes)+' genes:\n')
if len(taxaremoved)>0:
	for x in range(0,len(taxaremoved)):
		coveragefile.write(taxaremoved[x]+'\n')
coveragefile.write('\nTerminals per gene:\n')
for x in range(0,len(datasets)):
	nameholder = inputfilenames[x][:(len(inputfilenames[x])-len(fastaextension))]    # remove extension
	nameholder = nameholder.split('/') # get rid of path
	nameholder = nameholder[len(nameholder)-1]
	namefornexus = nameholder.translate(str.maketrans('-. /\\','_____')) # get rid of name parts that might trouble phylogenetics software
	coveragefile.write(namefornexus + '\t' + str(genecoverage[x]) + '\n')
coveragefile.write('\nGenes removed because less than '+str(mintaxa)+' terminals:\n')
if len(genesremoved)>0:
	for x in range(0,len(genesremoved)):
		nameholder = genesremoved[x][:(len(genesremoved[x])-len(fastaextension))]    # remove extension
		nameholder = nameholder.split('/') # get rid of path
		nameholder = nameholder[len(nameholder)-1]
		namefornexus = nameholder.translate(str.maketrans('-. /\\','_____')) # get rid of name parts that might trouble phylogenetics software
		coveragefile.write(namefornexus + '\n')
coveragefile.close()
