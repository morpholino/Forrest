from Bio import AlignIO
import argparse

parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-i', '--infile', help='Fasta/Phylip set to be trimmed', required=True)
parser.add_argument('-f', '--format', help='Alignment format', default="infer")


args = parser.parse_args()

if args.format == "infer":
	suffix = args.infile.split(".")[-1]
	if suffix in ["fasta", "fna", "faa", "fas", "fa"]:
		msaformat = "fasta"
	elif suffix in ["phylip", "phy"]:
		msaformat = "phylip-relaxed"
	elif suffix in ["nex", "nexus"]:
		msaformat = "nexus"
	else:
		quit("unrecognized MSA format")
elif args.format in ["fasta", "phylip", "clustal", "emboss", "nexus", "stockholm"]:
	msaformat = args.format
else:
	quit("unrecognized MSA format")


align = AlignIO.read(args.infile, msaformat)
filt = align[:, :0] #has to be two ranges!
#print(filt)
seqs, length = len(align), align.get_alignment_length()
#print(align)
print(seqs,length)
for i,col in enumerate(range(0,length)):
	if i % 1000 == 0:
		print("processed: {}".format(i))
	if align[:, col].count("-") == seqs:
		continue
	if align[:, col].count("X") < seqs:
		#print(align[:, col:col+1]) #just col won't work
		filt += align[:, col:col+1]
filt = [record for record in filt if record.seq.count("-") != len(record.seq) and len(record.seq) != 0]

#AlignIO.write(filt, args.infile.replace("." + suffix, "_noX." + suffix), msaformat)
outfile = args.infile.replace("." + suffix, "_noX." + suffix)
with open(outfile, "w") as result:
	for index, r in enumerate(filt):
		#get rid of the trailing newline character at the end of file:
		if index != len(filt) - 1:
			result.write(">{}\n{}\n".format(r.id, r.seq)) #r.description includes seqlen
		else:
			result.write(">{}\n{}".format(r.id, r.seq))

print("Alignment dimensions: {} seqs x{} pos".format(len(filt), len(filt[0].seq)))
