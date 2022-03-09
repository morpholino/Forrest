from Bio import SeqIO
import argparse

def _aggregate_union(data):
    sorted_data = sorted(data)
    result = sorted_data[0]

    for lower, upper in sorted_data[1:]:
        # If we ever find result[1] is None, we know it covers any
        # other possible values, so we can stop at that point.
        if result[1] is None:
            return None

        if lower > result[1]:
            print("Non-overlapping motif!")
            #keep the entire span or return the original motif coordinates
            #return result

        result = (result[0], max(upper, result[1]))

    return result


def _parse_gff(ingff, motif):
	ranges_d = {}
	print(motif)

	with open(ingff) as f:
		for l in f:
			if l.startswith("##FASTA"):
				break
			if l.startswith("#"):
				continue
			data = l.strip().split("\t")
			seqname = data[0]
			annotation = data[-1].lower()
			if any(x in annotation for x in motif):
				if seqname not in ranges_d:
					ranges_d[seqname] = []
				ranges_d[seqname].append([int(data[3]), int(data[4])])

	return ranges_d


def _parse_fasta(infasta):
	return {seq.name: seq.seq for seq in SeqIO.parse(infasta, "fasta")}


def trim_fasta(seq_d, ranges_d, out_trimmed, out_failed):
	# Sorts sequences into those having a conserved domain and those which do not
	# Seqs with conserved domains are trimmed
	t,f = 0,0
	with open(out_trimmed, "wt") as trimmed, \
		 open(out_failed, "wt") as failed:
		for seqname in seq_d:
			if seqname in ranges_d:
				motif = _aggregate_union(ranges_d[seqname])
				if motif:
					sequence = seq_d[seqname][motif[0]-1:motif[1]]
				else:
					print("ERROR: No range for", seqname, ranges_d[seqname])
				trimmed.write(">{}__{}-{}\n{}\n".format(seqname, motif[0], motif[1], sequence))
				t += 1
			else:
				failed.write(">{}\n{}\n".format(seqname, seq_d[seqname]))
				f += 1
	print("Passed: {}; failed: {}.".format(t,f))


def main():
	parser = argparse.ArgumentParser(description='How to use argparse')
	parser.add_argument('-f', '--fastain', help='Fasta/Phylip set to be trimmed', required=True)
	parser.add_argument('-g', '--gffin', help='GFF input file', required=True)
	parser.add_argument('-m', '--motif', help='Motif(s) to filter on', required=True)
	args = parser.parse_args()

	ingff = args.gffin #"HSP20_bact_v3.gff"
	infasta = args.fastain #"HSP20_bact_v3.fasta"
	suffix = infasta.split(".")[-1]
	out_trimmed = infasta.replace("." + suffix, "-trimmed." + suffix)
	out_failed = infasta.replace("." + suffix, "-failed." + suffix)
	motif = args.motif.lower().split(",") #"HSP20,HEAT SHOCK PROTEIN,SHSP"

	ranges_d = _parse_gff(ingff, motif)
	#print(ranges_d)
	seq_d = _parse_fasta(infasta)
	trim_fasta(seq_d, ranges_d, out_trimmed, out_failed)


if __name__ == '__main__':
	main()