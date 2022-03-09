import os

files = [x for x in os.listdir(".") if x.endswith("pdb.txt")]
codechange = {"ALA": "A", "ARG": "R", "ASN": "N",
"ASP": "D", "CYS": "C", "GLU": "E", 
"GLN": "Q", "GLY": "G", "HIS": "H",
"ILE": "I", "LEU": "L", "LYS": "K",
"MET": "M", "PHE": "F", "PRO": "P",
"SER": "S", "THR": "T", "TRP": "W",
"TYR": "Y", "VAL": "V", "HOH": ""}

with open("ribosomals.fasta", "w") as result:
	for f in files:
		name = f.replace("Suppl_Data_1_", "").replace(".pdb.txt", "")
		print(name)
		seq = 	""
		lastposition = "0"
		with open(f) as data:
			for l in data:
				position = l.split()[5]
				try:
					aa = codechange[l.split()[3]]
				except KeyError:
					print("unrecognized amino acid:", l.split()[3])
					aa = ""
				if position != lastposition:
					seq += aa
					lastposition = position

			print(seq)
			result.write(">{}\n{}\n".format(name, seq))

