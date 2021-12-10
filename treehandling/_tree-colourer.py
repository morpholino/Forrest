import os
import argparse
import sys
from Bio import SeqIO

print("This script colours branches of an input nexus tree based on presence/absence of motif names in leaf labels.")
print("usage: python tree-colourer.py tree.nex motifs.txt \n--------------------------------------------------------------------------")

outfile = sys.argv[1]
intree = open(sys.argv[1]).read()
motifs = open(sys.argv[2]).read().split()

basecolours_d = {'blue': '[&!color=#0000ff]', 
				 'brown': '[&!color=#996633]', 
				 'cyan': '[&!color=#00ffff]', 
				 'green': '[&!color=#00ff00]', 
				 'magenta': '[&!color=#ff00ff]', 
				 'orange': '[&!color=#ff8000]', 
				 'purple': '[&!color=#800080]', 
				 'red': '[&!color=#ff0000]' , 
				 'yellow': '[&!color=#ffff00] '}
				 
basecolours_l = list(basecolours_d.keys())
print(basecolours_l)
#colourcodes_l = ['[&!color=#0000ff]', '[&!color=#996633]', '[&!color=#00ffff]', '[&!color=#00ff00]', '[&!color=#ff00ff]', '[&!color=#ff8000]', '[&!color=#800080]', '[&!color=#ff0000]' , '[&!color=#ffff00]']
#black = ['-16777216', '000000']
motifcolours = {}
motifcolourcodes = {}
c = 0
for motif in motifs:
	c += 1
	currentcolour = basecolours_l[c]
	motifcolours[motif] = currentcolour
	motifcolourcodes[motif] = basecolours_d[currentcolour]

with open("motifcolours.txt", "w") as out:
	for key, value in motifcolours.items():
		print(key, ":\t", value)
		out.write(f"{key}:\t{value}")


#parse the nexus file
filestart = intree.split('\n')[:4]
taxa = ('\n'.join(intree.split('\n')[4:]).split(';')[0]).split('\n')
finaltaxa = list()
fileend = ';'.join('\n'.join(intree.split('\n')[4:]).split(';')[1:])

for line in taxa:
	for motif in motifs:
		if motif in line:
			line = line + motifcolourcodes[motif]
	finaltaxa.append(line)

# print(type(filestart))
# print(type(finaltaxa))
# print(type(fileend))

print("writing coloured tree as: coloured-{}".format(outfile))
#write results
with open('coloured-' + outfile, 'w') as out:
	for line in filestart:
		out.write(line +"\n")
	for line in finaltaxa:
		out.write(line +"\n")
	out.write(";"+fileend)