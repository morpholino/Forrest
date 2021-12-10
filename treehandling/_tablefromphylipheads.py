#for i in `ls *.linsi`; do trimal -in $i -out $i.trimal.phy -phylip -automated1; done
#head -n 1 *.phy > alignments.txt

infile = open("alignments.txt").read().split("\n")
outfile = open("alignments-table.tsv", "w")
outfile.write("gene\t\t\ttaxa\tlength\n")
for line in infile:
	if line.startswith("=="):
		line = line.split()[1]
		line = line.split(".")[0]
		outfile.write(line + "\t\t")
	elif len(line) > 0:
		line = line.split()
		outfile.write("\t\t".join(line) + "\n")
	else:
		pass