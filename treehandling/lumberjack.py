import os,re,argparse
from Bio import SeqIO

def clear_tag(nodename):
	partpattern = r'_([\.\d]+)%aligned'
	try:
		taghit = re.search(partpattern, nodename)
		tag = taghit.group()
		perc = taghit.group(1)
		nodename = nodename.replace(tag, "")
	except:
		perc = 101
		tag = ""
		#print(f"No pattern in {nodename}")
	return nodename

#set working directory
if os.path.isdir("/Users/morpholino/OwnCloud/"):
	home = "/Users/morpholino/OwnCloud/"
elif os.path.isdir("/Volumes/zoliq data/OwnCloud/"):
	home = "/Volumes/zoliq data/OwnCloud/"
else:
	print("Please set a homedir")

print("This script removes unwanted branches from a dataset based on an input nexus tree.")
print("Mark the branches to be removed in the input tree by a colour (using FigTree). Please use basic colours or format found in your nex file. ")
print("usage: python phylo-pruner.py -f eno.fasta -t testtree.nex [-p prefix_for_filtered -d working_directory -O -c all]")
print('for batches, move to dir of interest, then:\nfor i in *fasta; do j="${i%.fasta}" && python ~/OwnCloud/progs/PYTHON/treehandling/tree-reducer.py -f $i -t $j.treefile -p v2 -O; done\n')
#WIN:black = #-16777216, #000000; green = #-16737997, #009933


parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-f', '--fastain', help='Fasta/Phylip set to be trimmed', required=True)
parser.add_argument('-t', '--tree', help='Treefile for trimming', required=True)
parser.add_argument('-c', '--colour', help='Branch colour filter', default='all')
parser.add_argument('-p', '--prefix', help='Filtered dataset prefix', default='filtered')
parser.add_argument('-d', '--directory', help='Change working directory', default='.')
parser.add_argument('-O', '--no_omitted', help='Supress omitted', action='store_true')

args = parser.parse_args()
if args.directory == ".":
	#wd = home + "genomes/euglena longa/trees/todo"
	wd = os.getcwd()
	print(wd)
	os.chdir(wd)
else:
	os.chdir(args.directory)

if args.prefix == "-":
	prefix = "_"
else:
	prefix = args.prefix

print(args.fastain)
if args.fastain.startswith(prefix):
	quit("Target file already present? Quitting...\n{}".format(70*"-"))
if args.fastain.split(".")[-1] in ("fasta", "fas", "fst", "fa", "faa", "ali"):
	indataset = SeqIO.parse(args.fastain, 'fasta')
elif args.fastain.split(".")[-1] in ("phy", "phylip"):
	indataset = SeqIO.parse(args.fastain, 'phylip')
else:
	quit("file type not recognized - is it fasta/fas/fa/faa/fst or phy/phylip?")
intree = open(args.tree).read()
filtercolour = args.colour

#you can provide a homemade database renaming dictionary
taxarepl = {"actiCORYd": "Corynebacter diphteriae"}


basecolours = {'blue': '0000ff', 'brown': '996633', 'cyan': '00ffff', 'green': '00ff00', 'magenta': 'ff00ff', 'orange': 'ff8000', 'purple': '800080', 'red': 'ff0000' , 'yellow': 'ffff00', 'white': 'ffffff'}
black = ['-16777216', '000000']
if filtercolour in basecolours:
	filtercolour = basecolours[filtercolour]
elif filtercolour == 'all':
	print("any colour accepted")
else:
	print("unknown filter, setting to 'user-defined'. taxa with unrecognized colour codes will be retained")

#load fasta
seq_d = {}
badchars = ("|#=+,:;()'[]/") #also []/
for sequence in indataset:
	shortname = sequence.description.replace(" ","_")
	newname = []
	for c in shortname:
		if c in badchars:
			c = "_"
		newname.append(c)
	shortname = ''.join(newname)
	#print(shortname)
	seq_d[shortname] = sequence.seq
	if shortname.split("_")[0] in taxarepl:
		taxon = taxarepl[shortname.split("_")[0]]
		shortname = shortname.replace(shortname.split("_")[0], taxon.replace(" ", "_"))
		seq_d[shortname] = sequence.seq
	#print(">" + shortname + "\n" + seq_d[shortname])

print("done loading sequences")
#load taxa from tree
alltaxa = [] # a list of leaf labels with their color
pattern = r"\[&!color=.+"
skip = []
skippedc = 0
keptc = 0
treelines = intree.split('\n')[4:] #extract only leaf labels
for line in treelines:
	if line != ';':
		line = line.replace("'", "")
		line = line.replace("\t", "").replace(" ","_").replace("|","_")#.replace("@","_")
		if "%aligned" in line:
			line = clear_tag(line)
		#line = line.split("@")[0]
		#linecolour = line + colour
		#print(linecolour)
		alltaxa.append(line)
		if "[&!color" in line:
			colour = re.search(pattern, line).group()
			newcolour = line.split('[&!color=#')[1].replace("]", "")
			#print(taxon)
			if newcolour in black:
				print("black detected for %s, keeping this taxon" % (line))
				keptc += 1
			elif filtercolour == 'all':
				skip.append(line.split('[&!color=#')[0])
				skippedc += 1
			elif newcolour == filtercolour:
				skip.append(line.split('[&!color=#')[0])
				skippedc += 1
			else:
				print("unknown colour detected for %s, keeping this taxon" % (line))
				keptc += 1
			newtaxon = line.split('[&!color=#')[0]
			#print(newtaxon)
			alltaxa = [newtaxon if x==line else x for x in alltaxa] 
		else:
			#print(line)
			keptc += 1
			colour = ""

	else:
		break

errorfile = open("_key_errors.txt", "a")
print("done loading taxa")
if not args.no_omitted:
	print("omitted taxa listed in omitted-{}".format(args.fastain))
	#write omitted taxa
	with open('omitted-' + args.fastain, 'w') as f:
		for taxon in skip:
			nohighertaxon = taxon.split("@")[0]
			nospacetaxon = taxon.replace(" ","_")
			if seq_d.get(taxon) != None:
				f.write(">%s\n%s\n" % (taxon, seq_d[taxon]))
			elif seq_d.get(nohighertaxon) != None:
				f.write(">%s\n%s\n" % (taxon, seq_d[nohighertaxon]))
			elif seq_d.get(nospacetaxon) != None:
				f.write(">%s\n%s\n" % (taxon, seq_d[nospacetaxon]))
			else:
				print("!!!!!KEY ERROR for omitted", taxon)
				errorfile.write("omitted:\t{}\n".format(taxon))

print("writing filtered dataset...")
#write results
with open('{}-{}'.format(args.prefix, args.fastain), 'w') as out:
	for taxon in alltaxa:
		nohighertaxon = taxon.split("@")[0]
		nospacetaxon = taxon.replace(" ","_")
		if taxon not in skip:
			if seq_d.get(taxon) != None:
				out.write(">%s\n%s\n" % (taxon, seq_d[taxon]))
			elif seq_d.get(nohighertaxon) != None:
				out.write(">%s\n%s\n" % (taxon, seq_d[nohighertaxon]))
			elif seq_d.get(nospacetaxon) != None:
				out.write(">%s\n%s\n" % (taxon, seq_d[nospacetaxon]))
			else:
				print("!!!!!KEY ERROR for filtered", taxon)
				errorfile.write("filtered:\t{}\n".format(taxon))

print("WRITING DONE, \n\t{} taxa kept in {}-{},".format(keptc, args.prefix, args.fastain))
if not args.no_omitted:
	print("\t{} taxa omitted in omitted-{}".format(skippedc, args.fastain))
else:
	print("\t{} taxa omitted".format(skippedc))
print("Filtering finished!")
print("\n--------------------------------------------------------------------------")
