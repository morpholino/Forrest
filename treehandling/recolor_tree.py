import os,re,argparse
from Bio import SeqIO

# https://github.com/morpholino


def _makeset(labelfiles, labeltypes, len_threshold=200, trim_names=False):
	# Returns a set of seqIDs extracted from labelfiles
	# For gff files, returns any sequences having an annotation longer than len_threshold
	labelset = []

	labelfile = [x for x in labelfiles if any(l in x for l in labeltypes)]
	if len(labelfile) > 0:
		labelfile = labelfile[0]
	else:
		return labelset

	if os.path.exists(labelfile):
		if labelfile.split('.')[-1] in ('fasta', 'fas', 'fst', 'fa', 'faa', 'ali'):
			labelset = [x.name for x in SeqIO.parse(labelfile, 'fasta')]
			if trim_names:
				print("Trimming seqnames to remove __xx-xx")
				labelset = [x.split("__")[0] for x in labelset]
		elif labelfile.split('.')[-1] in ('tsv', 'txt'):
			with open(labelfile, 'rt') as file:
				labelset = {x.strip().split('\t')[0] for x in file.readlines()}
				labelset = list(labelset)
		elif labelfile.split('.')[-1] in ('gff', 'gff3'):
			labelset = set()
			with open(labelfile, 'rt') as file:
				for l in file:
					l = l.strip().split('\t')
					try:
						genelen = int(l[4]) - int(l[3])
						if genelen > len_threshold:
							labelset.add(l[0])
					except:
						print("not numbers", l[3], l[4])
			labelset = list(labelset)
		else:
			print("Unknown file extension, quitting!", labelfile)
	else:
		if labelfile == '':
			pass
		else:
			print('{} labels file {} not found!'.format(labeltype, labelfile))
	
	return labelset


def recolor_tree(intree, labelset, color):
	# Returns a nex tree as string, with tip labels marked with color
	treeidx = intree.find('begin trees;')
	#print(treeidx)
	if treeidx == -1:
		os.system("open {}".format(intree))
		quit("Tree not in nexus format!")
	taxablock = intree[:treeidx]
	treeblock = intree[treeidx:]

	basecolors = {'blue': '[&!color=#0000ff]', 
				 'brown': '[&!color=#996633]', 
				 'cyan': '[&!color=#00ffff]', 
				 'green': '[&!color=#00ff00]', 
				 'magenta': '[&!color=#ff00ff]', 
				 'orange': '[&!color=#ff8000]', 
				 'purple': '[&!color=#800080]', 
				 'red': '[&!color=#ff0000]' , 
				 'yellow': '[&!color=#ffff00]',
				 'black': '[&!color=#000000]'}

	if color not in basecolors:
		raise ValueError("Color not in base colors!")
	color_re = r'\[&!color=#\w+\]'
	for item in labelset:
		itemidx = taxablock.find(item)
		if itemidx == -1:
			print("Label not found!", item)
			continue

		if taxablock[itemidx+len(item)] == "'":
			item = item + "'"
		if taxablock[itemidx+len(item)] not in ("\n", "["):
			#but shorter IDs usually come first in nexus!
			print("imperfect match", item, itemidx)
			itemidx = taxablock.find(item, itemidx+1)
		newitem = "{}{}".format(item, basecolors[color])
		if itemidx > 0:
			length = len(item) + 17 #the length of a color annotation
			extendeditem = taxablock[itemidx:itemidx+length]
			if re.search(color_re, extendeditem):
				taxablock = taxablock.replace(extendeditem, newitem)
			else:
				#prevent shorter ID from replacing longer IDs!
				taxablock = taxablock.replace("{}\n".format(item), "{}\n".format(newitem))
		else:
			print("Label not found!", item)
			continue

	intree = taxablock + treeblock
	return intree


def main():
	#set working directory
	if os.path.isdir('/Users/morpholino/OwnCloud/'):
		home = '/Users/morpholino/OwnCloud/'
	elif os.path.isdir('/Volumes/zoliq data/OwnCloud/'):
		home = '/Volumes/zoliq data/OwnCloud/'
	else:
		print('Please set a homedir')

	#WIN:black = #-16777216, #000000; green = #-16737997, #009933
	print("Recolor trees")

	parser = argparse.ArgumentParser(description='How to use argparse')
	parser.add_argument('-t', '--tree', help='Treefile for trimming', required=True)
	parser.add_argument('-p', '--prefix', help='Filtered datasets prefix', default='-')
	parser.add_argument('-g', '--good_labels', help='File with good labels - will be marked green', default='')
	parser.add_argument('-b', '--bad_labels', help='File with bad labels - will be marked red', default='')
	parser.add_argument('-o', '--other_labels', help='File with unknown-identity labels - will be marked orange', default='')
	parser.add_argument('-d', '--directory', help='Change working directory', default='.')

	#parser.add_argument('-U', '--no_unknown', help='Supress unknowns', action='store_true')

	args = parser.parse_args()
	#Some more help:
	# This script marks branches from a nexus tree based on input lists.
	# Usage: python recolor_tree.py -t testtree.nex [-p prefix_for_filters -g good_labels -b bad_labels -o other_labels -d working_directory]
	# For batches, move to dir of interest, then:
	# for i in *dbhits.fasta; do j=${i%-dbhits.fasta}; python ~/OwnCloud/progs/PYTHON/treehandling/recolor_tree.py -p "$j"-dbhits -t fast-"$j".merged.tree; done
	# for i in *fasta; do j=${i%.fasta}; python ~/OwnCloud/progs/PYTHON/treehandling/recolor_tree.py -p "$j" -t "$j".tree; done

	if args.directory == '.':
		wd = os.getcwd()
		#wd = home + 'AndyLab/Phaeocystis/annotation/sulfur metabolism/MATOU'
		print(wd)
		os.chdir(wd)
	else:
		os.chdir(args.directory)

	if os.path.exists(args.tree):
		with open(args.tree, 'rt') as f:
			intree = f.read()
	else:
		quit('Tree file not found! Quitting...\n')

	if args.prefix == '-':
		if all([x == '' for x in [args.good_labels, args.bad_labels, args.other_labels]]):
			quit('No filter files given! Quitting...\n')
		else:
			labelfiles = [args.good_labels, args.bad_labels, args.other_labels]
	else:
		prefix = args.prefix
		#find all good, bad and unknown files with this prefix!
		labelfiles = [x for x in os.listdir('.') if x.startswith(prefix)]

	goodset = _makeset(labelfiles, ["good", "trimmed"], trim_names=True)
	badset = _makeset(labelfiles, ["bad", "failed"])
	badset = [x for x in badset if x not in goodset]
	unknownset = _makeset(labelfiles, ["unknown"], trim_names=True)

	if all(len(x) == 0 for x in [goodset, badset, unknownset]):
		quit("No filter files found among: {}! Quitting...\n".format(", ".join(labelfiles)))

	print("{}: good {}, bad {}, unknown {}".format(args.tree, len(goodset), len(badset), len(unknownset)))

	intree = recolor_tree(intree, goodset, "green")
	intree = recolor_tree(intree, badset, "red")
	intree = recolor_tree(intree, unknownset, "orange")

	print('writing recoloured tree as: recol-{}\n'.format(args.tree))
	#write results
	with open('recol-{}'.format(args.tree), 'w') as out:
		out.write(intree)


if __name__ == "__main__":
	main()