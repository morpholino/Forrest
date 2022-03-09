import csv
from Bio import SeqIO

with open('renaming_MMETSP.txt', 'r') as f:
	reader = csv.reader(f, delimiter='\t')
	mmetsp2species = {r[0]: r[1].replace(" ", "-") for r in reader}

c = 0
infasta = SeqIO.parse("MMETSP_pep.fa", "fasta")
with open("MMETSP_new.fa", "w") as result:
	for seq in infasta:
		c += 1
		if c % 100000 == 0:
			print(c, " processed")
		if seq.name.startswith("MMETSP"):
			mmetsp = seq.name[:10]
			try:
				taxon = mmetsp2species[mmetsp]
				result.write(">{}_{}\n{}\n".format(taxon, seq.description, seq.seq))
			except KeyError:
				print(seq.name, "taxon not retrieved")
				result.write(">{}\n{}\n".format(seq.description, seq.seq))
		else:
			result.write(">{}\n{}\n".format(seq.description, seq.seq))