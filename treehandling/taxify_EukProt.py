import argparse,os,gzip
import seaborn as sns
from Bio import SeqIO,Entrez
from ete3 import NCBITaxa
#http://etetoolkit.org/docs/2.3/tutorial/tutorial_ncbitaxonomy.html
ncbi = NCBITaxa()
Entrez.email = 'zoltan.fussy@natur.cuni.cz'
Entrez.api_key = "ed51bca6c73792ecc692af11ff762b72a008"
#update at times:
#ncbi.update_taxonomy_database()


def _get_cmap(n, name='Spectral'): #hsv for very divergent data?
	import matplotlib.pyplot as plt
	import matplotlib.colors as mcolors
	'''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
	RGB color; the keyword argument name must be a standard mpl colormap name.'''
	colormap = plt.cm.get_cmap(name, n)
	rgbcolors = []
	for i in range(colormap.N):
		rgb = colormap(i)[:3] # will return rgba, we take only first 3 so we get rgb
		rgbcolors.append(mcolors.rgb2hex(rgb)) #to hex
		#rgb = [int(256*x) for x in rgb]
		#rgbcolors.append(rgb) # to 16-bit rgb
	return rgbcolors
 

def _colormap(categories, palette):
	cmap = _get_cmap(len(categories), name=palette)
	colors = {}
	for i, X in enumerate(categories):
		colors[X] = cmap[i] #hopefully this is still recognized as a color
	return colors


def _find_homedir(directory):
	if os.path.isdir("/Users/morpholino/OwnCloud/"):
		home = "/Users/morpholino/OwnCloud/"
	elif os.path.isdir("/Volumes/zoliq data/OwnCloud/"):
		home = "/Volumes/zoliq data/OwnCloud/"
	else:
		print("Please set a homedir")

	if directory == ".":
		print("changing to default directory")
		default = "progs/forrest/treehandling/"
		return home + default
	else:
		return home + directory


def _dictionaries():
	#try using taxids to remove the need to translate the lineage into rank
	global manual_taxids
	manual_taxids = {2461416: ['Bacteria', 'Terrabacteria group', 'Cyanobacteria/Melainabacteria group', 'Cyanobacteria', 'unclassified Cyanobacteria'],
					 3152: ['Eukaryota', 'Viridiplantae', 'Chlorophyta', 'Chloropicophyceae', 'Chloropicales'],
					 2853422: ['Eukaryota', 'Malawimonadida']}

	global eukprot_supergroups
	eukprot_supergroups = {'Alveolata': 'Alveolata', 
						   'Amoebozoa': 'Amoebozoa', 
						   'Ancoracysta': 'basal-euks', 'Ancyromonadida': 'basal-euks', 'Collodictyonidae': 'basal-euks', 'Hemimastigophora': 'basal-euks', 
						   'Malawimonadida': 'basal-euks', 'Rotosphaerida': 'basal-euks', 
						   'Mantamonas': 'basal-euks', 'Rigifilida': 'basal-euks', 
						   'Telonemia': 'basal-euks', 'Telonemida': 'basal-euks',
						   'Apusomonadida': 'basal-opis', 'Apusozoa': 'basal-opis', 'Breviatea': 'basal-opis', 'Filasterea': 'basal-opis', 'Choanoflagellata': 'basal-opis', 
						   'Ichthyosporea': 'basal-opis', 'Pluriformea': 'basal-opis', 'Opisthokonta': 'basal-opis',
						   'Chloroplastida': 'Chloroplastida', 'Viridiplantae': 'Chloroplastida', 
						   'Cryptophyceae': 'Cryptista', 'Kathablepharidacea': 'Cryptista', 'Palpitomonas': 'Cryptista', 
						   'Euglenozoa': 'Discoba', 'Heterolobosea': 'Discoba', 'Jakobida': 'Discoba', 'Tsukubamonadida': 'Discoba',
						   'Fornicata': 'Metamonada', 'Parabasalia': 'Metamonada', 'Preaxostyla': 'Metamonada', 
						   'Fungi': 'Fungi', 
						   'Glaucophyta': 'Glaucophyta', 'Glaucocystophyceae': 'Glaucophyta',
						   'Phaeocystales': 'Haptista-PHC', 'Isochrysidales': 'Haptista-ISO', 'Prymnesiales': 'Haptista-PRY', 'Coccolithales': 'Haptista-COC',
						   'Haptophyta': 'Haptista', 'Centroplasthelida': 'Haptista', 
						   'Metazoa': 'Metazoa', 
						   'Rhodophyta': 'Rhodophyta', 'Picozoa': 'Rhodophyta', 'Rhodelphis': 'Rhodophyta', 
						   'Rhizaria': 'Rhizaria', 
						   'Stramenopiles': 'Stramenopiles',
						   }

	global eukprot_groups
	eukprot_groups = {'Alveolata': 'Alveolata', 'Dinoflagellata': 'Dinoflagellata', 'Apicomplexa': 'Apicomplexa', 'Ciliophora': 'Ciliophora',
					   'Amoebozoa': 'Amoebozoa', 
					   'Ancoracysta': 'basal-euks', 'Ancyromonadida': 'basal-euks', 'Collodictyonidae': 'basal-euks', 'Hemimastigophora': 'basal-euks', 
					   'Malawimonadida': 'basal-euks', 'Rotosphaerida': 'basal-euks', 
					   'Mantamonas': 'basal-euks', 'Rigifilida': 'basal-euks', 
					   'Telonemia': 'Telonemida', 'Telonemida': 'Telonemida',
					   'Apusomonadida': 'Apusozoa', 'Apusozoa': 'Apusozoa', 'Breviatea': 'Breviatea', 'Filasterea': 'Holozoa', 'Choanoflagellata': 'Choanozoa', 
					   'Ichthyosporea': 'Holozoa', 'Pluriformea': 'Holozoa', 'Opisthokonta': 'Opisthokonta',
					   'Opisthosporidia': 'Holomycota', 'Nucleariida': 'Holomycota', 'Fonticulida': 'Holomycota',
					   'Chlorophyta': 'Viridiplantae', 'Prasinodermophyta': 'Viridiplantae','Streptophyta': 'Streptophyta', 
					   'Cryptophyceae': 'Cryptista', 'Kathablepharidacea': 'Cryptista', 'Palpitomonas': 'Cryptista', 
					   'Euglenozoa': 'Euglenozoa', 'Heterolobosea': 'Heterolobosea', 'Jakobida': 'Jakobida', 'Tsukubamonadida': 'Tsukubamonadida',
					   'Fornicata': 'Fornicata', 'Parabasalia': 'Parabasalia', 'Preaxostyla': 'Preaxostyla', 
					   'Fungi': 'Fungi', 
					   'Glaucophyta': 'Glaucophyta', 'Glaucocystophyceae': 'Glaucophyta',
					   #'Phaeocystales': 'Haptista-PHC', 'Isochrysidales': 'Haptista-ISO', 'Prymnesiales': 'Haptista-PRY', 'Coccolithales': 'Haptista-COC',
					   'Haptophyta': 'Haptista', 'Centroplasthelida': 'Haptista', 
					   'Metazoa': 'Metazoa', 
					   'Rhodophyta': 'Rhodophyta', 'Picozoa': 'Rhodophyta', 'Rhodelphis': 'Rhodophyta', 
					   'Rhizaria': 'Rhizaria', 
					   'Stramenopiles': 'Stramenopiles',
					   }

	global bacterial_supergroups
	bacterial_supergroups = {'Acidobacteria': 'other-bct', 'Calditrichaeota': 'other-bct', 'Nitrospinae/Tectomicrobia group': 'other-bct', 
							 'Spirochaetes': 'other-bct', 'Synergistetes': 'other-bct', 'Coprothermobacterota': 'other-bct', 
							 'unclassified Bacteria': 'other-bct', 'Bacteria incertae sedis': 'other-bct', 'Nitrospirae': 'other-bct', 
							 'Thermodesulfobacteria': 'other-bct', 'Deferribacteres': 'other-bct', 'Thermotogae': 'other-bct', 'Aquificae': 'other-bct', 
							 'PVC group': 'PVC', 'FCB group': 'FCB-Bacteroidetes', 
							 # subdivisions for Proteobacteria
							 'Alphaproteobacteria': 'Alphaproteobacteria', 'Betaproteobacteria': 'Betaproteobacteria', 'Gammaproteobacteria': 'Gammaproteobacteria',
							 'delta/epsilon subdivisions': 'Deltaproteobacteria',
							 # subdivisions for Terrabacteria
							 'Actinobacteria': 'Actinobacteria', 'Cyanobacteria/Melainabacteria group': 'Cyanobacteria', 'Firmicutes': 'Firmicutes',
							 }

	global viral_supergroups
	viral_supergroups = {'Alphasatellitidae': 'virus', 'Riboviria': 'virus', 'Monodnaviria': 'virus', 'Naldaviricetes': 'virus',
						 #both NCLDV and Polinton viruses belong to Varidnaviria:
						 'Varidnaviria': 'Varidnaviria',
						 'Duplodnaviria': 'Duplodnaviria'
						 }


def _read_table(infile):
	outdict = {}
	if infile.endswith(".gz"):
		with gzip.open(infile, 'rt') as f:
			for l in f:
				l = l.strip().split("\t")
				if len(l) > 1:
					if l[1] in ("Taxid", "N/A"):
						continue
					outdict[l[0]] = int(l[1])
	else:
		with open(infile, 'rt') as f:
			for l in f:
				l = l.strip().split("\t")
				if len(l) > 1:
					if l[1] in ("Taxid", "N/A"):
						continue
					outdict[l[0]] = int(l[1])
	return outdict


def _translate_taxids(taxids_d):
	_dictionaries()
	eukaryotes = eukprot_groups
	clade_d = {}
	# 1 = root
	# 131567 = cellular organisms
	# 12908 = unclassified
	# 28384 = other (synthetic, transposons, etc.)
	undefined = {1, 131567, 12908, 28384}
	for key, taxid in taxids_d.items():
		supergroup = ''
		if taxid in undefined:
			continue

		try:
			if taxid in manual_taxids:
				rank = manual_taxids[taxid]
			else:
				lineage = ncbi.get_lineage(taxid)[2:]
				names = ncbi.get_taxid_translator(lineage)
				rank = [names[taxid] for taxid in lineage]
			supergroup = [x for x in rank if x in eukaryotes]
		except ValueError:
			print(f'problem with {key}:{taxid}, please try updating taxa.sqlite')
			with open('errors.log', 'a') as errorfile:
				errorfile.write(f'error retrieving taxid for {key}:{taxid}\n')
			continue

		try:
			# if sequence is bacterial
			if rank[0] == 'Bacteria':
				if len(rank) > 1:
					if rank[1] in ['Terrabacteria group', 'Proteobacteria']:
						if len(rank) > 2:
							supergroup = bacterial_supergroups.get(rank[2], 'other-bct')
						else:
							#print('Insufficient data: {}'.format(rank))
							supergroup = 'other-bct'
					else:
						supergroup = bacterial_supergroups.get(rank[1], 'other-bct')
				else:
					supergroup = 'other-bct'
				#supergroups.append('_'.join(rank[1:2]))

			# if archaeal
			elif rank[0] == 'Archaea':
				supergroup = 'Archaea'

			# if viral
			elif any(x in rank[0] for x in {'viria', 'satellitidae', 'viricetes'}):
				if rank[0] == 'Varidnaviria':
					if rank[2] == 'Nucleocytoviricota':
						supergroup = 'Varidnaviria-NCLDV'
					elif rank[2] == 'Preplasmiviricota':
						supergroup = 'Varidnaviria-Poli'
					else:
						print('Unknown varidnavirus {}'.format(rank))
				else:
					supergroup = viral_supergroups.get(rank[0], 'virus')
				#supergroups.append(rank[0]+'_'+rank[2])
			elif 'viruses' in rank[0]:
				if 'Pleurochrysis sp. endemic' in rank[1]:
					#both Pleurochrysis sp. endemic NCLDV and Polinton virophage
					supergroup = 'Varidnaviria-NCLDV'
				elif 'Pleurochrysis sp. Polinton' in rank[1]:
					supergroup = 'Varidnaviria-Poli'
				else:
					supergroup = 'virus'
				#supergroups.append(rank[0])

			# else eukaryotic
			elif len(supergroup) == 1:
				supergroup = eukaryotes.get(supergroup[0], 'Eukaryota')
			elif len(supergroup) > 1:
				#use the lowest major rank
				supergroup = eukaryotes.get(supergroup[-1], 'Eukaryota')
			#some unwanted situations
			elif rank == ['Eukaryota']:
				supergroup = 'Eukaryota'
			elif supergroup == []:
				print('No usable taxon data: {} => {}'.format(taxid, rank))
				continue
			else:
				print('Other error: {} => {}'.format(taxid, rank))
				continue
			clade_d[key] = supergroup
		except IndexError:
			print("Could not parse taxid", taxid, rank)

	return clade_d


def write_itol_ranges(infile, prefix, taxonomyfile):
	"""

	Input:		TSV with seqIDs in the first column or FASTA file

	Returns:	Writes an ITOL color_range file

	"""
	#this is for an experimental function -- automated color assignment from a spectral colormap
	"""taxa =  ['Rhodophyta',
						'empty',
						'Haptista',
						'empty',
						'Cryptista',
						'empty',
						'Alveolata',
						'Dinoflagellata',
						'Apicomplexa',
						'Ciliphora',
						'empty',
						'Rhizaria',
						'empty',
						'Stramenopiles',
						'empty',
						'Parabasalia',
						'Preaxostyla',
						'Fornicata',
						'empty',
						'Euglenozoa',
						'empty',
						'Streptophyta',
						'Viridiplantae',
						'empty',
						'Glaucophyta',
						'empty',
						'Amoebozoa',
						'Holomycota',
						'Fungi',
						'Metazoa',
						'Holozoa',
						'Choanozoa',
						'Breviatea',
						'Apusozoa',
						'Opisthokonta']
			
				bw =	[
						'basal-euks',
						'Jakobida',
						'Heterolobosea',
						'Telonemida',
						'Tsukubamonadida',
						'Eukaryota']"""

	color_map_ranges = {
						'Alveolata': '#FCAD67',
						'Amoebozoa': '#7AD3E5',
						'Apicomplexa': '#FC8822',
						'Apusozoa': '#A9C6C6',
						'basal-euks': '#B6BBBB',
						'Breviatea': '#A9BFC8',
						'Viridiplantae': '#47C662',
						'Choanozoa': ' #6A9EB4',
						'Ciliphora': '#D4BB8D',
						'Cryptista': '#95507A',
						'Dinoflagellata': '#E57D1C',
						'Euglenozoa': '#5ADABD',
						'Eukaryota': '#919191',
						'Fornicata': '#19A39A',
						'Fungi': '#68A6CF',
						'Glaucophyta': '#d2f6b6',
						'Haptista': '#90669E',
						'Heterolobosea': '#28867E',
						'Holomycota': '#87C0E6',
						'Holozoa': '#718EB0',
						'Jakobida': '#1D8C8A',
						'Metazoa': '#4F719A', 
						'Opisthokonta': '#7ABCE5',
						'Parabasalia': '#DAD3B1',
						'Preaxostyla': '#B7AF8B',
						'Rhizaria': '#8E7564',
						'Rhodophyta': '#BC3241',
						'Stramenopiles': '#B07A57',
						'Streptophyta': '#0BA42C',
						'Telonemida': '#D3BE65',
						'Tsukubamonadida': '#8CEAD5'
						}

	header = "TREE_COLORS\nSEPARATOR TAB\nDATA\n"
			#First 3 fields define the node, type and color
			#Possible node types are:
			#'range': defines a colored range (colored background for labels/clade => for taxonomic labelling)
			#'clade': defines color/style for all branches in a clade
			#'branch': defines color/style for a single branch
			#'label': defines font color/style for the leaf label
			#'label_background': defines the leaf label background color

	if os.path.isfile(taxonomyfile) == True:
		pass
	elif os.path.isfile("{}_taxids.tsv".format(prefix)) == True:
		taxonomyfile = "{}_taxids.tsv".format(prefix)
	else:
		print("Prepare *taxids.tsv first using the --writetaxids option!")
	taxid_d = _read_table(taxonomyfile)
	print("Imported taxids:", len(taxid_d.keys()))

	print("Loading list of seqIDs...")
	if infile.split(".")[-1] in ("fasta", "faa", "fas"):
		#script can process fasta as infile; note that FASTA have modified seqnames!
		seqids = {x.name for x in SeqIO.parse(infile, "fasta")}
	else:
		with open(infile, "rt") as f:
			seqids = {x.split("\t")[0] for x in f}

	#this should contain clade info for EukProt codes
	clade_dict = _translate_taxids(taxid_d)
	with open("{}_taxids.txt".format(prefix), "wt") as taxonfile:
		taxonfile.write(header)
		for seqid in seqids:
			clade = clade_dict.get(seqid.split("_")[0], "")
			taxonfile.write("{}\trange\t{}\t{}\n".format(seqid, color_map_ranges.get(clade, "#DCDCDC"), clade))
	
	#experimental automatic color assignment
	"""color_map_automatic = _colormap(taxa, palette="Spectral")
				for t in taxa:
					if t == "empty":
						continue
					hex = color_map_automatic[t]
					print(f"{t} {hex}")"""


def write_taxids(infile, prefix, taxonomyfile):
	"""Compares a subset of the MATOU database that contains unigenes with the given Pfam

	Input:		TSV with seqIDs in the first column or FASTA file

	Returns:	Writes a file

	"""
	print("Loading list of seqIDs...")
	if infile.split(".")[-1] in ("fasta", "faa", "fas"):
		#script can process fasta as infile; note that FASTA have modified seqnames!
		seqids = {x.name for x in SeqIO.parse(infile, "fasta")}
	else:
		with open(infile, "rt") as f:
			#assuming seq name is in first column
			seqids = {x.split("\t")[0] for x in f}

	if taxonomyfile.endswith(".gz"):
		with gzip.open(taxonomyfile, "rt") as f,\
			open("{}_taxids.tsv".format(prefix), "wt") as taxonfile:
			for l in f:
				seqid = l.split("\t")[0]
				if seqid in seqids:
					taxonfile.write("{}\t{}\n".format(seqid, l.split("\t")[1]))
	else:
		with open(taxonomyfile, "rt") as f,\
			open("{}_taxids.tsv".format(prefix), "wt") as taxonfile:
			for l in f:
				seqid = l.split("\t")[0]
				if seqid in seqids:
					taxonfile.write("{}\t{}\n".format(seqid, l.split("\t")[1]))



def taxify_file(infile, prefix, taxonomyfile, translate=True):
	"""Makes dictionaries for translating seqIDs to taxids or clades.

	Input:		TSV with seqIDs in the first column and taxid in the second
				File to modify

	Returns:	Writes a modified file

	"""
	if os.path.isfile(taxonomyfile) == True:
		pass
	elif os.path.isfile("{}_taxids.tsv".format(prefix)) == True:
		taxonomyfile = "{}_taxids.tsv".format(prefix)
	else:
		print("Prepare *taxids.tsv first using the --writetaxids option!")
	taxid_d = _read_table(taxonomyfile)
	print("Imported taxids:", len(taxid_d.keys()))
	suffix = infile.split(".")[-1]
	#Taxids are in raw form - only numbers, therefore in tree and fasta files 
	# we need to extend their names.
	if translate:
		outfile = infile.replace(".{}".format(suffix), ".clades.{}".format(suffix))
		clade_d = _translate_taxids(taxid_d)
		with open(infile, 'rt') as f:
			file = f.read()
			for seqid in clade_d:
				file = file.replace("{}".format(seqid), "{}.{}".format(clade_d[seqid], seqid))
		with open(outfile, 'wt') as result:
			result.write(file)
		#print("Translated taxids:", clade_d)
	else:
		outfile = infile.replace(".{}".format(suffix), ".taxids.{}".format(suffix))
		with open(infile, 'rt') as f:
			file = f.read()
			for seqid in taxid_d:
				file = file.replace("{}".format(seqid), "{}.{}".format(taxid_d[seqid], seqid))
		with open(outfile, 'wt') as result:
			result.write(file)



def main():
	parser = argparse.ArgumentParser(description='How to use argparse')
	#use a Pfam ID or a list/fasta of MATOU sequences
	group = parser.add_mutually_exclusive_group()
	parser.add_argument('-t', '--taxify', help='Apply taxids in this file (fasta or tree)', default='')
	group.add_argument('-i', '--infile', help='Infile to be used to subset taxids (fasta or list)', default='')
	group.add_argument('-p', '--prefix', help='Infile prefix', default='')
	parser.add_argument('--taxonomy', help='TSV file containing custom taxonomy', default='')
	parser.add_argument('--translate', help='Translate taxids into clades', action='store_true')
	parser.add_argument('--writetaxids', help='Write a file of taxids', action='store_true')
	parser.add_argument('-w', '--writeitol', help='Write ITOL color_ranges', action='store_true')

	args = parser.parse_args()
	if args.infile != "":
		infile = args.infile
		prefix = infile.split(".")[0]
	elif args.prefix != "":
		infile = args.prefix + ".list"
		prefix = args.prefix
	else:
		#neither prefix nor prefix list was given
		prefix = "tmp"

	#if os.path.isfile(infile) == False:
	#	os.system('zgrep "{0}" prefix.tsv.gz > {0}.list'.format(prefix))

	if args.taxonomy == "":
		taxonomyfile = "{}/EukProt_included_data_sets.v02.taxids.tsv".format(_find_homedir("genomes"))
	else:
		taxonomyfile = args.taxonomy

	if args.writetaxids:
		write_taxids(infile, prefix, taxonomyfile)

	if args.writeitol:
		write_itol_ranges(infile, prefix, taxonomyfile)

	if args.taxify != "":
		taxify_file(args.taxify, prefix, taxonomyfile, args.translate)


if __name__ == '__main__':
	main()
