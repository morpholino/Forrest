import os,argparse
from Bio import SeqIO
import gzip

#os.chdir("/Users/zoliq/ownCloud/progs/PYTHON/mydata")
"""
prefixes = [x for x in os.listdir("./alignments") if x.endswith(".ali")]
#prefixes = ["PF00134_rp75", "PF02984_rp75"]
for prefix in prefixes:
	#os.system("hmmbuild {0}.hmm {0}.ali".format(prefix))
	os.system("hmmsearch -o {0}.result.txt --notextw --tblout {0}.table.txt --cpu 2 {0}.hmm hmmdb.fa".format(prefix, prefix))
"""
def _parse_seqdb(fastadb, allseqids):
	seq_d = {}
	if fastadb.endswith("gz"):
		with gzip.open(fastadb, "rt") as handle:
			for seq in SeqIO.parse(handle, "fasta"):
				if seq.name in allseqids:
					seq_d[seq.name] = str(seq.seq).replace("*", "")			
	else:
		for seq in SeqIO.parse(fastadb, "fasta"):
			if seq.name in allseqids:
				seq_d[seq.name] = str(seq.seq).replace("*", "")
	return seq_d


def _parse_hmmsearch_file(files, seq_d, ignore_threshold):
	for file in files:
		multidomain = set()
		written = set() #move multidomain, written into the for loop to have domain-specific files
		name = file.replace(".table.txt", "-dbhits.fasta")
		with open(file) as infile, open(name, "w") as result:
			for line in infile:
				if "Domain annotation for each sequence" in line:
					break
				if "inclusion threshold" in line:
					if ignore_threshold:
						continue
					else:
						break
				line = line.split()
				if len(line) < 9:
					continue
				if  line[0] != "#":
					try:
						evalue = float(line[0])
						seqname = line[8]
						#print(seqname)
						try:
							if "::" in seqname:
								#process in-house transcriptome seqIDs to be shorter
								writeseqname = "{}__{}_{}".format(seqname.split("::")[0], seqname.split("::")[-2], seqname.split("::")[-1])
							else:
								writeseqname = seqname
						except IndexError:
							writeseqname = seqname
							print("Could not parse seqname", writeseqname)
						if seqname in written:				
							pass
						elif evalue < 1: #which is always...
							result.write(">{}\n{}\n".format(writeseqname, seq_d[seqname]))
							written.add(seqname)
							multidomain.add(seqname)
						elif seqname in multidomain: #several domains found in the same seq
							result.write(">{}\n{}\n".format(writeseqname, seq_d[seqname]))
							written.add(seqname)
						else:
							multidomain.add(seqname)
					except ValueError:
						pass
						#this is not a score table line


def hmmparser():
	parser = argparse.ArgumentParser(description='How to use argparse')
	parser.add_argument('-d', '--database', help='input database', default="EukProt_v2_renamed.faa")
	parser.add_argument('-s', '--hmmsearch', help='HMMsearch output', default="batch")
	parser.add_argument('-i', '--ignore_threshold', help='Ignore inclusion threshold in hmmer output', action="store_true")
	args = parser.parse_args()

	if os.path.isfile(args.database):
		fastadb = args.database
	else:
		quit("{} database not found!".format(args.database))
	
	#read in seqs from DB fasta only if seqname in results:
	if args.hmmsearch == "batch":
		files = [x for x in os.listdir(".") if x.endswith("table.txt")]
	else:
		files = args.hmmsearch.split(",")
	files = [x for x in files if os.path.isfile(x)]
	if files == []:
		quit("HMMsearch input files not found!")
	
	allseqids = set()
	for file in files:
		data = [x for x in open(file).readlines() if x.startswith(">>")]
		for item in data:
			allseqids.add(item.split()[1]) #this should be seqID even if there is a description
	if allseqids == {}:
		quit("No hits in hmm input(s)!")

	seq_d = _parse_seqdb(fastadb, allseqids)

	#filter results and export fasta
	_parse_hmmsearch_file(files, seq_d, args.ignore_threshold)


if __name__ == "__main__":
	hmmparser()