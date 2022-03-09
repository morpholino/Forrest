from Bio import SeqIO
import argparse
import gzip
import time
import re
from ete3 import Tree


def fasta_rename(infile, outfile, outtable, keyfile, keycolumn, add_prefix):
	print("Fasta file identified as input, processing {} > {}".format(infile, outfile))
	seq_dict = {}
	for seq in SeqIO.parse(infile, "fasta"):
		seq_dict[seq.name] = {"seq": seq.seq, "newname": seq.name}

	found = 0
	if keyfile.endswith(".gz"): #gzip.open is much slower and does not work?
		start = time.time()
		with gzip.open(keyfile, "rt") as f:
			for l in f:
				l = l.strip().split("\t")
				if l[0] in seq_dict:
					found += 1
					if add_prefix:
						newname = "{}__{}".format("_".join(l[0].split("_")[1:3]), l[keycolumn])
					else:
						newname = l[keycolumn]
					seq_dict[l[0]]["newname"] = newname
		end = time.time()
	else:
		start = time.time()
		with open(keyfile, "rt") as f:
			for l in f:
				l = l.strip().split("\t")
				if l[0] in seq_dict:
					found += 1
					if add_prefix:
						newname = "{}__{}".format("_".join(l[0].split("_")[1:3]), l[keycolumn])
					else:
						newname = l[keycolumn]
					seq_dict[l[0]]["newname"] = newname
		end = time.time()

	print("Time elapsed: {}".format(end - start))
	print("Found {} keys in dictionary (out of {} sequences)".format(found, len([x for x in seq_dict.keys() if x.startswith("EP0")])))	

	with open(outfile, "wt") as result:
		for seq in seq_dict:
			result.write(">{}\n{}\n".format(seq_dict[seq]["newname"], seq_dict[seq]["seq"]))

	with open(outtable, "wt") as result:
		for key, item in seq_dict.items():
			if not key.startswith("EP"):
				continue
			result.write("{}\t{}\n".format(key, item["newname"]))


def newick_rename(infile, outfile, outtable, keyfile, keycolumn, add_prefix):
	print("Tree file identified as input, processing {} > {}".format(infile, outfile))
	branch_dict = {}
	support_re = r'(\[&label=([\d\.]+)\])'
	with open(infile) as f:
		treedata = f.read()
	if treedata.startswith("#NEXUS"):
		print("Nexus format detected, switching to nexus_rename()!")
		nexus_rename(infile, outfile, outtable, keyfile, keycolumn, add_prefix)
		quit("Finished!")
	for hit in re.findall(support_re, treedata):
		treedata = treedata.replace(hit[0], hit[1])
	#print(treedata)
	t = Tree(treedata)
	print("File open")

	all_leaves = t.get_tree_root().get_leaves()
	new_leaves = {x.name: x.name for x in all_leaves}
	print("{}, {} from EukProt".format(len(all_leaves), len([x for x in all_leaves if x.name.startswith("EP0")])))
	found = 0
	if keyfile.endswith(".gz"): #gzip.open is somewhat slower
		start = time.time()
		with gzip.open(keyfile, "rt") as f:
			for l in f:
				l = l.strip().split("\t")
				if l[0] in new_leaves:
					found += 1
					if add_prefix:
						newname = "{}__{}".format("_".join(l[0].split("_")[1:3]), l[keycolumn])
					else:
						newname = l[keycolumn]
					new_leaves[l[0]] = newname
		end = time.time()
	else:
		start = time.time()
		with open(keyfile, "rt") as f:
			for l in f:
				l = l.strip().split("\t")
				if l[0] in new_leaves:
					found += 1
					if add_prefix:
						newname = "{}__{}".format("_".join(l[0].split("_")[1:3]), l[keycolumn])
					else:
						newname = l[keycolumn]
					new_leaves[l[0]] = newname
		end = time.time()
	print("Time elapsed: {}".format(end - start))
	print("Found {} keys in dictionary (out of {} sequences)".format(found, len([x for x in all_leaves if x.name.startswith("EP0")])))	

	#print(new_leaves)
	for node in t.traverse():
		if node.is_leaf():
			node.name = new_leaves[node.name]
	t.write(format=0, outfile=outfile)

	with open(outtable, "wt") as result:
		for key, item in new_leaves.items():
			if not key.startswith("EP"):
				continue
			result.write("{}\t{}\n".format(key, item))


def nexus_rename(infile, outfile, outtable, keyfile, keycolumn, add_prefix):
	#1 read file into string
	print("Tree file identified as input, processing {} > {}".format(infile, outfile))
	with open(infile) as f:
		filedata = f.read()
	if not filedata.startswith("#NEXUS"):
		quit("Not nexus format!")

	#2 extract tree block and leaves using ete3 as for nexus
	treeblock_re = r'begin trees;\n\ttree tree_1 = \[&R\] (.*)\nend;'
	#treeblock_re = r'begin trees;(.*)end;'
	support_re = r'(\[&label=([\d\.]+)\])'
	try:
		treedata = re.search(treeblock_re, filedata).group(1)
	except AttributeError:
		treedata = filedata.split("begin trees;\n")[1].split("\nend;")[0]
		if len(treedata.split("\n")) > 1:
			quit("More than one tree found, unsupported file arrangement!")
	for hit in re.findall(support_re, treedata):
		treedata = treedata.replace(hit[0], hit[1])
	#print(treedata)
	t = Tree(treedata)
	print("File open")

	branch_dict = {}
	badchars = ["::", ","]
	all_leaves = t.get_tree_root().get_leaves()
	new_leaves = {x.name: x.name for x in all_leaves}
	print("{}, {} from EukProt".format(len(all_leaves), len([x for x in all_leaves if x.name.startswith("EP0")])))
	found = 0
	if keyfile.endswith(".gz"): #gzip.open is somewhat slower
		start = time.time()
		with gzip.open(keyfile, "rt") as f:
			for l in f:
				l = l.strip().split("\t")
				if l[0] in new_leaves:
					found += 1
					if add_prefix:
						newname = "{}__{}".format("_".join(l[0].split("_")[1:3]), l[keycolumn])
					else:
						newname = l[keycolumn]
					if any(x in newname for x in badchars):
						for x in badchars:
							newname = newname.replace(x, "_")
					new_leaves[l[0]] = newname
		end = time.time()
	else:
		start = time.time()
		with open(keyfile, "rt") as f:
			for l in f:
				l = l.strip().split("\t")
				if l[0] in new_leaves:
					found += 1
					if add_prefix:
						newname = "{}__{}".format("_".join(l[0].split("_")[1:3]), l[keycolumn])
					else:
						newname = l[keycolumn]
					if any(x in newname for x in badchars):
						for x in badchars:
							newname = newname.replace(x, "_")
					new_leaves[l[0]] = newname
		end = time.time()
	print("Time elapsed: {}".format(end - start))
	print("Found {} keys in dictionary (out of {} sequences)".format(found, len([x for x in all_leaves if x.name.startswith("EP0")])))	

	#3 sort new_leaves dictionary according to string length - longer seq ID's will be replaced first!
	leaves_list = [x for x in new_leaves.keys()]
	leaves_list.sort(key=len, reverse=True)

	#4 loop through new_leaves and replace original file data
	for l in leaves_list:
		filedata = filedata.replace(l, new_leaves[l])

	with open(outfile, "wt") as result:
		result.write(filedata)

	with open(outtable, "wt") as result:
		for key, item in new_leaves.items():
			if not key.startswith("EP"):
				continue
			result.write("{}\t{}\n".format(key, item))


def main():
	parser = argparse.ArgumentParser(description='How to use argparse')
	parser.add_argument('-i', '--infile', help='File to process', required=True)
	parser.add_argument('-o', '--outfile', help='Output filename', default="")
	parser.add_argument('-k', '--keyfile', help='Treefile for trimming', default="EukProt_v2_rename.key")
	parser.add_argument('--column', help='Keyfile column', default=1)
	parser.add_argument('--add_prefix', help='Add species prefix', action='store_true')

	args = parser.parse_args()

	infile = args.infile
	suffix = infile.split(".")[-1]
	if args.outfile == "":
		outfile = infile.replace(suffix, "renamed.{}".format(suffix))
		outtable = infile.replace(suffix, "rename.tsv")
	else:
		outfile = args.outfile
		outtable = outfile.replace(suffix, "tsv")
	keyfile = args.keyfile
	keycolumn = args.column

	if suffix in ("fasta", "fas", "fna", "faa", "fa"):
		fasta_rename(infile, outfile, outtable, keyfile, keycolumn, args.add_prefix)

	elif suffix in ("tree", "treefile", "nwk"):
		newick_rename(infile, outfile, outtable, keyfile, keycolumn, args.add_prefix)

	elif suffix in ("nexus", "nex", "nxs"):
		nexus_rename(infile, outfile, outtable, keyfile, keycolumn, args.add_prefix)
		#maybe a different method for nexus trees?


if __name__ == '__main__':
	main()
