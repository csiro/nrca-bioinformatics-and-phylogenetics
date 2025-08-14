import getopt, string, sys, glob, os

outputfile = './structure_summary.txt'
pathtogenes = './'

#parse parameters

try:
	opts, args = getopt.getopt(sys.argv[1:],"hi:s:")
except getopt.GetoptError:
	print('structure_summary.py -i <inputdirectory> -s <nameofanalysis>')
	print('Name of analysis must not have any extensions, e.g., if the fastStructure output is blah.3.meanQ for K = 3, then the file name must be blah.')
	print('More details: structure_summary.py -h')
	sys.exit(2)
for opt, arg in opts:
	if opt == '-h':
		print('Script for summarising results of a fastStructure analysis across several values of K.')
		print('\nUse: structure_summary.py -i <inputdirectory> -s <nameofanalysis>')
		print('\nParameters')
		print('h: Help, i.e. get this information')
		print('i: Path to folder containing files to be evaluated. Default is current working directory.')
		#print('o: Output file name for summary. Default is structure_summary.txt in current working directory')
		print('s: Analysis name. This must not have any extensions, e.g., if the fastStructure output is blah.3.meanQ for K = 3, then the file name must be blah.')
		sys.exit()
	elif opt in ("-i"):
		pathtogenes = arg
		if pathtogenes[len(pathtogenes)-1] != "/":
			pathtogenes = pathtogenes + "/"
	elif opt in ("-s"):
		analysisname = arg

inputfilenames = glob.glob(pathtogenes + analysisname + ".*.log")
ks = [int(os.path.basename(k).split('.')[1]) for k in inputfilenames]
mls = ['']*(max(ks)-min(ks)+1)
ksord = ks
ksord.sort()
kstxt = [str(k) for k in ksord]

for inputfile in inputfilenames:
	with open(inputfile, 'r') as f:
		data = f.read()
	f.close()
	data = data.split('\n')
	thisk = int(os.path.basename(inputfile).split('.')[1])
	for thisline in data:
		if thisline[0:len("Marginal Likelihood = ")] == "Marginal Likelihood = ":
			thisml = thisline[len("Marginal Likelihood = "):]
	mls[thisk-min(ks)] =  thisml

#write output file

outfile = open(pathtogenes+analysisname+'.summary', "w")
for k in range(0,len(ks)):
	outfile.write(kstxt[k]+"\t"+mls[k]+"\n")
outfile.close()

