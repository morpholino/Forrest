import re
import sys

orgnpattern = r"OS=(.*)OX="

with open(sys.argv[1]) as f:
	result = f.read()
	#print(result)
suffix = sys.argv[1].split(".")[-1]
outfile = sys.argv[1].replace(suffix, "uni."+suffix)

rawresult = result.replace("tr|", "").replace("sp|", "")
resultstring = ""
for line in rawresult.split(">")[1:]:
	try:
		orgn = re.search(orgnpattern, line).group(1)
		orgn = "_".join(orgn.split()[:2])
		acc_entry = line.split()[0].replace("|", "@")
		annot = " ".join(line.split(" ")[1:]).split("OS=")[0]
		seq = "".join(line.split("\n")[1:])
		resultstring += f">{orgn}_{acc_entry} {annot}\n{seq}\n"
	except IndexError:
		pass

with open(outfile, "w") as out:
	out.write(resultstring)