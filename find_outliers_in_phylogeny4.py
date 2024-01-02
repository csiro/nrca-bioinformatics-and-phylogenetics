from Bio import Phylo
import sys
import csv

# only parameter is tree file name
treefilename = sys.argv[1]

# read tree file

myphylogeny = Phylo.read(treefilename, "newick")

# tree already has family names on it, so no need to have a separate taxonomy file

# figure out which family each terminal belongs to

treeterminalsaffiliation =[]
treeterminals = myphylogeny.get_terminals()
#print(treeterminals)

treeterminalsname = [i.name for i in treeterminals]
for x in range(0,len(treeterminalsname)):
	treeterminalsaffiliation.append((treeterminalsname[x].split('_'))[0])
families = list(set(treeterminalsaffiliation))
#print(families)

# get all internal nodes from phylogeny

internalnodes = myphylogeny.get_nonterminals()

# for each family, we need to find out what the node is with the maximum of
# ( descendants of node that belong to family / number of family members ) * ( descendants of node that belong to family / total number of descendants of node )
# seems easiest to build up lists, one for each family,
# by appending a number for each node as we walk through them

monophylyscores = [[] for n in range(len(families))]

# first, count members in each family

familymembers = [0 for n in range(len(families))]

for x in range(0,len(treeterminalsaffiliation)):
	familymembers[families.index(treeterminalsaffiliation[x])] += 1

# loop through internal nodes and score for each family

for node in range(0,len(internalnodes)):
	cladeterminals = internalnodes[node].get_terminals()
	familymembersinclade = [0 for n in range(len(families))]
	for trm in range(0,len(cladeterminals)):
		familymembersinclade[families.index(treeterminalsaffiliation[treeterminalsname.index(cladeterminals[trm].name)])] += 1
	for fam in range(0,len(families)):
		newscore = (float(familymembersinclade[fam]) / float(familymembers[fam])) * (float(familymembersinclade[fam]) / float(len(cladeterminals)))
		monophylyscores[fam].append(newscore)

# find out most representative nodes for each family

bestnodes = []
for fam in range(0,len(families)):
	bestnodes.append(monophylyscores[fam].index(max(monophylyscores[fam])))

# now go through each family and identify its members that are not descendants of its most represenative node

outliers=[]
for fam in range(0,len(families)):
	if familymembers[fam] > 0:
		cladeterminals = internalnodes[bestnodes[fam]].get_terminals()
		cladeterminalnames = set([i.name for i in cladeterminals])
		familymembersoverall = set([treeterminalsname[i] for i in range(len(treeterminalsname)) if treeterminalsaffiliation[i]==families[fam]])
		newoutliers = list(familymembersoverall.difference(cladeterminalnames))
		if len(newoutliers) > 0:
			for x in range(0,len(newoutliers)):
				outliers.append(newoutliers[x])

listfile = open(treefilename+'.outliers.txt','w')
if len(outliers) > 0:
	for x in range(0,len(outliers)):
		listfile.write(outliers[x]+'\n')
else:
	listfile.write('no outliers found\n')
listfile.close()
