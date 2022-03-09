from Bio import AlignIO
import sys

infile = sys.argv[1]
outfile = sys.argv[2]

alignment = AlignIO.read(infile, "stockholm")

with open(outfile, "wt") as handle:
	AlignIO.write(alignment, handle, "fasta")