#!/usr/bin/python3
import os
import re
from Bio import Entrez,SeqIO
from ete3 import NCBITaxa
#http://etetoolkit.org/docs/2.3/tutorial/tutorial_ncbitaxonomy.html
ncbi = NCBITaxa()

Entrez.email = 'zoltan.fussy@gmail.com'
#Entrez.api_key = os.environ["API_KEY"]
#print(os.environ)

if os.path.isdir("/Users/morpholino/OwnCloud/"):
	home = "/Users/morpholino/OwnCloud/"
elif os.path.isdir("/Volumes/zoliq data/OwnCloud/"):
	home = "/Volumes/zoliq data/OwnCloud/"
else:
	print("Please set a homedir")

#wd = "VoboraLab/data/publikace/140929 euglena/Volume3"
#wd = "AndyLab/Phaeocystis/annotation/metE"
#wd = "AndyLab/mTOR"
wd = "HamplLab/retortamonas/mito_protein_queries"
#wd = "DolezalLab/SecY/"
os.chdir(home + wd)

#will need these regexes
taxonpattern = r'\[(.+)\]'
genbankpattern = r'[A-Z]{1,3}(_)?\d+(\.\d)?'
badchars = ("|@+,:;()'")

#collect files
fastas = [x for x in os.listdir(".") if x.endswith("_out.fasta")] #NOTE: default output suffix of the diamondparse script
fastas.sort()

#main
for file in fastas:
	print("=============\n\n\nNow parsing file", file)
	fasta = SeqIO.parse(file, "fasta")
	with open(file.replace("_out.fasta", ".faa"), 'w') as out: #"_v1.fasta"
		for f in fasta:
			fname = f.name.strip()
			fdesc = f.description
			prot_idregex = re.search(genbankpattern, fname) #prot_idregex = without group(0)
			taxonregex = re.search(taxonpattern, fdesc)
			#this is for cases where taxon pattern is in the description
			if prot_idregex and taxonregex:
				prot_id = prot_idregex.group(0)
				taxon = "_".join(taxonregex.group(1).split()[:2])
				desc = f'{taxon}_{prot_id} {fdesc.replace(fname, "").split(" [")[0]}'
				seq = f.seq
			
			#this is for cases that have only gi|xxx|gb|yyyy| 
			elif fname.startswith("gi|"):
				gid = fname.split("|")[1]
				#print("GenBank ID found, yay!")
				print("Fetching data", gid)
				try:
					prot = Entrez.efetch(db='protein', id=gid, rettype='fasta', retmode='text')
					prot_record = SeqIO.read(prot, 'fasta')
					desc = prot_record.description
					taxon = re.search(taxonpattern, desc).group(1)
					if taxon:
						taxon = "_".join(taxon.split()[:2])
						desc = f'{taxon}_{desc.split(" [")[0]}'
					seq = prot_record.seq
					print("success", desc, seq[:20], "...")
				except:
					#bad request
					print("Bad request")
					desc = fdesc #fname
					seq = f.seq

			#this is for cases when XP_xxxx type AC is in seqname
			elif prot_idregex:
				prot_id = prot_idregex.group(0)
				#print("GenBank ID found, yay!", prot_id)
			#https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly
				if fname.startswith(prot_id):
					print("Fetching data", prot_id)
					try:
						prot = Entrez.efetch(db='protein', id=prot_id, rettype='fasta', retmode='text')
						prot_record = SeqIO.read(prot, 'fasta')
						desc = prot_record.description
						taxon = re.search(taxonpattern, desc).group(1)
						if taxon:
							taxon = "_".join(taxon.split()[:2])
							desc = f'{taxon}_{desc.split(" [")[0]}'
						seq = prot_record.seq
						print("success", desc, seq[:20], "...")
					except:
						#bad request
						print("Bad request")
						desc = f.name
						seq = f.seq
				else:
					desc = f.description
					seq = f.seq
			else:
				desc = f.description
				seq = f.seq
			out.write(f">{desc}\n{seq}\n")
			#out.write('{}\t{}__{}\n'.format(prot_id, orgn, tax))


quit("\n\nFinished processing, yay!")
#this is rich format:
for prot_id in ids:
	#print(prot_id)
	prot = Entrez.efetch(db='protein', id=prot_id, rettype='gb', retmode='text')
	prot_record = SeqIO.read(prot, 'genbank')
	print(prot_record)
	#tax = str(prot_record.annotations['taxonomy'][::-1]).replace('\'', '')
	orgn = prot_record.annotations['organism']
	name2taxid = ncbi.get_name_translator(orgn)
	print(prot_id, orgn, name2taxid)
	#out.write('{}\t{}__{}\n'.format(prot_id, orgn, tax))