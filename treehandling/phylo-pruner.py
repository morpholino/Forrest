import os,re,argparse
from Bio import SeqIO

def clear_tag(nodename):
	partpattern = r'_([\.\d]+)%aligned'
	try:
		taghit = re.search(partpattern, nodename)
		tag = taghit.group()
		perc = taghit.group(1)
		nodename = nodename.replace(tag, "")
	except:
		perc = 101
		tag = ""
		#print(f"No pattern in {nodename}")
	return nodename

#set working directory
if os.path.isdir("/Users/morpholino/OwnCloud/"):
	home = "/Users/morpholino/OwnCloud/"
elif os.path.isdir("/Volumes/zoliq data/OwnCloud/"):
	home = "/Volumes/zoliq data/OwnCloud/"
else:
	print("Please set a homedir")

print("This script removes unwanted branches from a dataset based on an input nexus tree.")
print("Mark the branches to be removed in the input tree by a colour (using FigTree). Please use basic colours or format found in your nex file. ")
print("usage: python phylo-pruner.py -f eno.fasta -t testtree.nex [-p prefix_for_filtered -d working_directory -O -c all]")
print('for batches, move to dir of interest, then:\nfor i in *fasta; do j="${i%.fasta}" && python ~/OwnCloud/progs/PYTHON/treehandling/tree-reducer.py -f $i -t $j.treefile -p v2 -O; done\n')
#WIN:black = #-16777216, #000000; green = #-16737997, #009933


parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-f', '--fastain', help='Fasta/Phylip set to be trimmed', required=True)
parser.add_argument('-t', '--tree', help='Treefile for trimming', required=True)
parser.add_argument('-c', '--colour', help='Branch colour filter', default='all')
parser.add_argument('-p', '--prefix', help='Filtered dataset prefix', default='filtered')
parser.add_argument('-d', '--directory', help='Change working directory', default='.')
parser.add_argument('-O', '--no_omitted', help='Supress omitted', action='store_true')

args = parser.parse_args()
if args.directory == ".":
	#wd = home + "genomes/euglena longa/trees/todo"
	wd = os.getcwd()
	print(wd)
	os.chdir(wd)
else:
	os.chdir(args.directory)

if args.prefix == "-":
	prefix = "_"
else:
	prefix = args.prefix

print(args.fastain)
if args.fastain.startswith(prefix):
	quit("Target file already present? Quitting...\n{}".format(70*"-"))
if args.fastain.split(".")[-1] in ("fasta", "fas", "fst", "fa", "faa", "ali"):
	indataset = SeqIO.parse(args.fastain, 'fasta')
elif args.fastain.split(".")[-1] in ("phy", "phylip"):
	indataset = SeqIO.parse(args.fastain, 'phylip')
else:
	quit("file type not recognized - is it fasta/fas/fa/faa/fst or phy/phylip?")
intree = open(args.tree).read()
filtercolour = args.colour

#zoltan's homemade database renaming dictionary
taxarepl9 = {"actiCORYd": "Corynebacter diphteriae", "actiMYCOt": "Mycobacterium tuberculosis", 
"actiSTREc": "Streptomyces coelicolor", "alfaAZOSs": "Azospirillum sp. B506", 
"alfaMAGNm": "Magnetospirillum magneticum", "alfaNITRh": "Nitrobacter hamburgensis", 
"alfaRHODs": "Rhodobacter sphaeroides", "alfaRICKc": "Rickettsia conori", 
"amoACANc": "Acanthamoeba castellanii str neff", "amoCAVEf": "Cavenderia fasciculata", 
"amoDICTd": "Dictyostelium discoideum", "amoDICTp": "Dictyostelium purpureum", 
"amoENTAh": "Entamoeba histolytica", "amoFILAn": "Filamoeba nolandi", "amoPLANf": "Planoprotostelium fungivorum",
"amoPOLYp": "Polysphondylium pallidum pn500", "amoSEXAN": "Sexangularia sp. CB-2014", 
"amoSTYGA": "Stygamoeba regulata", "amoTIEGl": "Tieghemostelium lacteum", "amoVANEL": "Vannella", 
"apiBABEb": "Babesia bovis", "apiCHROv": "Chromera velia", "apiCRYPm": "Cryptosporidium muris", 
"apiEIMEt": "Eimeria tenella", "apiGREGn": "Gregarina niphandrodes", "apiNEOSc": "Neospora caninum", 
"apiPLASb": "Plasmodium berghei", "apiPLASg": "Plasmodium gallinaceum",
"apiTHEIe": "Theileria equi", "apiTOXOg": "Toxoplasma gondii", "apiVITbr": "Vitrella brassicaformis", 
"apiVORp": "Voromonas pontica", "apuTHECt": "Thecamonas trahens",
"archPICRt": "Picrophilus torridus", "archPYROa": "Pyrobaculum aerophilum", 
"archSULFt": "Sulfolobus tokodaii", "archTHERv": "Thermoplasma volcanium", "bcidBACTf": "Bacteroides fragilis", 
"bcidFLAVc": "Flavobacterium columnare", "bcidPORPg": "Porphyromonas gingivalis", 
"bcidPREVr": "Prevotella ruminicola", "betaBURKc": "Burkholderia cenocepacia", "betaCUPRn": "Cupriavidus necator", 
"betaRALSs": "Ralstonia solanacearum", "betaVERMe": "Verminephrobacter eiseniae", 
"chlaAMORC": "Amorphochlora amoebiformis", "chlaBIGna": "Bigelowiella natans", "chlaCHLOR": "Chlorarachnion reptans", 
"chlaLOTam": "Lotharella amoebiformis", "chlaLOTgl": "Lotharella globosa", "chlaLOTgZ": "Lotharella globosa", 
"chlaLOTHs": "Lotharella sp. CCMP622", "chlaLOToc": "Lotharella oceanica", "chlaPARTg": "Partenskyella glossopodia", 
"chryCHATs": "Chattonella subsalsa CCMP2191", "chryDINOB": "Dinobryon sp UTEXLB2267", 
"chryHETak": "Heterosigma akashiwo CCMP3107", "chryOCHRO": "Ochromonas sp. CCMP1393", 
"chryVAUCH": "Vaucheria litorea CCMP2940", 
"cilCLIMv": "Climacostomum virens", "cilICHTm": "Ichthyophthirius multifiliis","cilMESOD": "Mesodinium pulex", 
"cilLITOp": "Litonotus pictus","cilOXYTt": "Oxytricha trifallax",
"cilPARAt": "Paramecium tetraurelia", "cilPLATm": "Platyophrya macrostoma", "cilPROTa": "Protocruzia adherens", 
"cilPSEUp": "Pseudocohnilembus persalinus", "cilSTENc": "Stentor coeruleus", "cilTETRt": "Tetrahymena thermophila", 
"crypCHROO": "Chroomonas cf. mesostigmatica", "crypCRYpa": "Cryptomonas paramecium", "crypCRYsp": "Cryptophyceae sp. CCMP2293",
"crypCRYPp": "Cryptomonas paramecium", "crypGONIp": "Goniomonas pacifica", "crypGUILt": "Guillardia theta", 
"crypGUIth": "Guillardia theta", "crypPALPb": "Palpitomonas bilix", "crypRHODl": "Rhodomonas lens", 
"cyanANABv": "Anabaena variabilis", "cyanCROCw": "Crocosphaera watsonii", "cyanCYANp": "Cyanothece sp. PCC 7425", 
"cyanLYNGp": "Lyngbya sp. PCC 8106", "cyanNODUs": "Nodularia spumigena", "cyanPROCm": "Prochlorococcus marinus", 
"cyanSYNEs": "Synechococcus sp. PCC 7335", "cyanTHERe": "Thermosynechococcus elongatus", 
"dinALEXa": "Alexandrium andersonii CCMP2222", "dinALEXc": "Alexandrium catenella OF101", 
"dinALEXt": "Alexandrium tamarense", "dinAMPHc": "Amphidinium carterae", "dinCERAf": "Ceratium fusus PA161109", 
"dinCRYPc": "Crypthecodinium cohnii WH Provasly-Seligo", "dinCRYPZ": "Crypthecodinium cohnii", 
"dinDURIb": "Durinskia baltica CSIRO CS-38", "dinDURIZ": "Durinskia baltica", "dinDINac": "Dinophysis acuminata",
"dinDYNac": "Dinophysis acuminata",
"dinGLENf": "Glenodinium foliaceum CCAP 1116-3", "dinGYMNc": "Gymnodinium catenatum", 
"dinHETro": "Heterocapsa rotundata SCCAP K-0483", "dinHETtr": "Heterocapsa triquestra CCMP 448", 
"dinKARb": "Karenia brevis SP1", "dinKARL": "Karlodinium veneficum", "dinKARmi": "Karlodinium micrum CCMP2283", 
"dinKARZ": "Karenia brevis", "dinKRYPf": "Kryptoperidinium foliaceum", "dinLINGp": "Lingulodinium polyedra CCMP 1738", 
"dinNOCTs": "Noctiluca scintillans", "dinOXYma": "Oxyrrhis marina", "dinPERKm": "Perkinsus marinus", 
"dinPRORm": "Prorocentrum minimum", "dinPYROb": "Pyrodinium bahamense", "dinSCRIP": "Scrippsiella hangoei", 
"dinSYMBs": "Symbiodinium sp.", "dinSYMBa": "Symbiodinium microadriaticum",
"eugANGOd": "Angomonas deanei", "eugBODOs": "Bodo saltans", 
"eugEUGgr": "Euglena gracilis", "eugEUGlo": "Euglena longa", "eugEUTc": "Eutreptiella gymnastica-like CCMP1594", 
"eugEUTn": "Eutreptiella gymnastica NIES-381", "eugEUTgn": "Eutreptiella gymnastica NIES-381", "eugLEISm": "Leishmania major", 
"eugLEPTp": "Leptomonas pyrrhocoris", "eugLEPTs": "Leptomonas seymouri",
"eugNEOBd": "Neobodo designis", "eugPARco": "Paratrypanosoma confusum", "eugPHYTO": "Phytomonas sp em1", 
"eugTRYPb": "Trypanosoma brucei", 
"excADUpa": "Aduncisulcus paluster", "excBLATn": "Blattamonas nauphoetae", "excCARPm": "Carpediemonas membranifera", 
"excCHIca": "Chilomastix caulleri", "excCHIcu": "Chilomastix cuspidata", "excDYSNb": "Dysnectes brevis", 
"excERGcy": "Ergobibamus cyprinoides", "excGIARi": "Giardia intestinalis P15", "excHISTm": "Histomonas meleagridis", 
"excIOTAs": "Iotanema sp.", "excKIPFb": "Kipferlia bialata", "excMONOe": "Monocercomonoides exilis", "excNAEgr": "Naegleria gruberi", 
"excPERco": "Percolomonas cosmopolitus ATCC50343", "excPTRIp": "Paratrimastix pyriformis", 
"excSPIRs": "Spironucleus salmonicida ATCC50377","excTREPs": "Trepomonas sp. PC1", "excTRIMm": "Trimastix marina", 
"extTTRIf": "Tritrichomonas foetus",
"firmBACIa": "Bacillus anthracis", "firmLISTm": "Listeria monocytogenes", "firmSTAPa": "Staphyllococcus aureus", 
"funASPEf": "Aspergillus fumigatus", "funCRYPn": "Cryptococcus neoformans", "funDEBAh": "Debaryomyces hansenii cbs767", 
"funLACCb": "Laccaria bicolor", "funNEURc": "Neurospora crassa", "funPUCCg": "Puccinia graminis", 
"funSCHIp": "Schizosaccharomyces pombe", 
"gamaSHEWb": "Shewanella baltica", "gamaVIBRc": "Vibrio cholerae", "gamaYERSp": "Yersinia pestis", 
"glauGLOEw": "Gloeochaete wittrockiana", "glauCYApa": "Cyanophora paradoxa", 
"glauCYPTg": "Cyanoptyche gloeocystis SAG4.97", "glauGLOwi": "Gloeochaete wittrockiana SAG46.84", 
"grnASTEs": "Asterochloris sp.", "grnAUXEp": "Auxenochlorella protothecoides",
"grnBATHp": "Bathycoccus prasinos", "grnCHLAl": "Chlamydomonas leiostraca", "grnCHLAl": "Chlamydomonas reinhardtii",
"grnCHLOs": "Chlorella sp.", "grnCHROz": "Chromochloris zofingiensis",
"grnCOCCs": "Coccomyxa subellipsoidea", "grnDUNAs": "Dunaliella salina", 
"grnDUNAt": "Dunaliella tertiolecta", 
"grnGONIp": "Gonium pectorale", "grnHELIs": "Helicosporidium sp.", 
"grnMICco": "Micromonas commoda", "grnMICpu": "Micromonas pusilla CCMP1545", "grnMICRZ": "Micromonas pusilla", 
"grnMONOn": "Monoraphidium neglectum", "grnNEPHp": "Nephroselmis pyriformis", "grnOSTRl": "Ostreococcus lucimarinus",
"grnOSTRm": "Ostreococcus mediterraneus", "grnOSTRt": "Ostreococcus tauri", "grnPICCL": "Picochlorum sp.", 
"grnPICCY": "Picocystis salinarum", "grnPOLYp": "Polytomella parva", "grnPOLYZ": "Polytomella parva", 
"grnPRASc": "Prasinococcus capsulatus", "grnPYRam": "Pyramimonas amylifera CCMP720", 
"grnPYRpa": "Pyramimonas parkeae CCMP726", "grnPYRpZ": "Pyramimonas parkeae", "grnTETRa": "Tetraselmis astigmatica", 
"grnTETRs": "Tetraselmis striata", "grnVOLVc": "Volvox carteri f. nagariensis", "hapCALCl": "Calcidiscus leptoporus", 
"hapCHRYt": "Chrysochromulina tobin", 
"hapEXANg": "Exanthemachrysis gayraliae", "hapIMANT": "Imantonia sp.", "hapISOCHg": "Isochrysis galbana", 
"hapPAVLs": "Pavlovales sp. CCMP2435", "hapPHAan": "Phaeocystis antarctica", "hapPHAgl": "Phaeocystis globosa",
"hapPHEOC": "Phaeocystis sp.", "hapPLEUc": "Pleurochrysis carterae", "hapPRYMp": "Prymnesium parvum", 
"hapSCYPH": "Scyphosphaera apsteinii", "hapEMILh": "Emiliania huxleyi", "haptEMIZ": "Emiliania huxleyi", 
"haptPLEUc": "Pleurochrysis carterae", "haptPRYMp": "Prymnesium parvum Texoma1", 
"hetPERCO": "Percolomonas cosmopolitus", "kytAMBOt": "Amborella trichopoda", "kytARABt": "Arabidopsis thaliana", 
"grnARABl": "Arabidopsis lyrata", "grnCHARb": "Chara braunii", 
"kytGLYCm": "Glycine max", "kytHORDv": "Hordeum vulgare", "kytORYZs": "Oryza sativa Japonica", 
"kytPHYSp": "Physcomitrella patens", "grnPHYSp": "Physcomitrella patens", "grnPOPtr": "Populus trichocarpa",
"kytSELAm": "Selaginella moellendorffii", "grnSELAm": "Selaginella moellendorffii", "grnSORGb": "Sorghum bicolor",
"kytSOlyc": "Solanum lycopersicum", "kytZEAma": "Zea mays", 
"metANOLc": "Anolis carolinensis", "metCAENe": "Caenorhabditis elegans", 
"metDAPHp": "Daphnia pulex", "metHELOr": "Helobdella robusta", "metLOAlo": "Loa loa",
"metLOTTg": "Lottia gigantea", 
"metNEMAv": "Nematostella vectensis", "metSCHIm": "Schistosoma mansoni", "metSTROp": "Strongylocentrotus purpuratus", 
"metTRICa": "Trichoplax adhaerens", "opiACANT": "Acanthoeca", "opiCAPSo": "Capsaspora owczarzaki",
"opiFONTa": "Fonticula alba",
"opiMONbr": "Monosiga brevicollis", "opiMONOb": "Monosiga brevicollis", 
"opiSALPr": "Salpingoeca rosetta", 
"redAHNFf": "Ahnfeltiopsis flabelliformis", "redBETAp": "Betaphycus philippinensis", 
"redCHONc": "Chondrus crispus", "redCOMPS": "Compsopogon caeruleus", "redCYANm": "Cyanidioschyzon merolae", 
"redERYTa": "Erythrolobus australicus", "redERYTm": "Erythrolobus madagascarensis", 
"redERYTm": "Erythrolobus madagascarensis", "redGALDs": "Galdieria sulphuraria", 
"redGRACv": "Gracilaria vermiculophylla", "redGRATt": "Grateloupia turuturu", "redHETSp": "Heterosiphonia pulchra",
"redMADAe": "Madagascaria erythrocladioides", "redMADAe": "Madagascaria erythrocladioides", 
"redPORae": "Porphyridium aerugineum SAG 1380-2", "redPORIp": "Porphyridium purpureum", 
"redPORum": "Porphyra umbilicalis", "redPORpu": "Porphyra purpurea", "redRHOma": "Rhodella maculata CCMP736", 
"redRHOSO": "Rhodosorus marinus", "redTIMSP": "Timspurckia oligopyrenoides", 
"rhiAMMOs": "Ammonia sp.", "rhiAMMOZ": "Ammonia sp.", 
"rhiELPHm": "Elphidium margitaceum", "rhiMINch": "Minchinia chitonis", "rhiPARTg": "Partenskyella glossopodia RCC365", 
"rhiPAULc": "Paulinella chromatophora", "rhiPLASb": "Plasmodiophora brassicae", "rhiRETIf": "Reticulomyxa filosa",
"rhiROSAL": "Rosalina", "rhiSORIT": "Sorites sp.", 
"strACHLh": "Achlya hypogyna", "strALBUl": "Albugo laibachii", "strAPHAi": "Aphanomyces invadans", 
"strAPLAN": "Aplanochytrium", "strAURAl": "Aurantiochytrium limacinum", "strAURAZ": "Aurantiochytrium limacinum", 
"strAUREN": "Aureococcus anophagefferens", "strAUREE": "Aureococcus anophagefferens", 
"strAUREa": "Aureococcus anophagefferens", "strAUREl": "Aureoumbra lagunensis", 
"strBICOS": "Bicosoecida sp. CB-2014", "strBLASh": "Blastocystis hominis", "strBOLID": "Bolidomonas sp.", 
"strBOLIp": "Bolidomonas pacifica", "strCAFro": "Cafeteria roenbergensis", "strCAFsp": "Cafeteria sp.", 
"strCAFca": "Cafeteria str Caron", "strCHATs": "Chattonella subsalsa", "strCHRRH": "Chrysoreinhardia sp.", 
"strCYCLc": "Cyclotella cryptica",
"strCYLIN": "Cylindrotheca closterium", "strDETON": "Detonula confervacea", "strDICTY": "Dictyocha speculum", 
"strDINOB": "Dinobryon sp.", "strDITYL": "Ditylum brightwellii", "strECTsi": "Ectocarpus siliculosus", 
"strEXTUs": "Extubocellulus spinifer", "strFRAGc": "Fragilariopsis cylindrus", 
"strFRAGk": "Fragilariopsis kerguelensis", "strFISTs": "Fistulifera solaris",
"strHETak": "Heterosigma akashiwo", 
"strHONDf": "Hondaea fermentalgiana", "strHYALa": "Hyaloperonospora arabidopsidis",
"strMALLO": "Mallomonas sp.", "strNANNg": "Nannochloropsis gaditana", "strNANNo": "Nannochloropsis oceanica", 
"strNITZs": "Nitzschia sp.", 
"strOCHRO": "Ochromonas sp.", "strOCHR2": "Ochromonadaceae sp. CCMP2298",
"strODONa": "Odontella aurita", "strPARAb": "Paraphysomonas bandaiensis", 
"strPELAG": "Pelagomonas calceolata", "strPELAs": "Pelagophyceae CCMP2097", "strPHAtr": "Phaeodactylum tricornutum", 
"strPHEOM": "Phaeomonas parva", 
"strPHYin": "Phytophtora infestans", "strPHYca": "Phytophthora capsici", "strPHYci": "Phytophthora cinnamomi",
"strPHYpa": "Phytophthora parasitica", "strPHYra": "Phytophtora ramorum", "strPHYso": "Phytophtora sojae", 
"strPINGU": "Pinguiococcus pyrenoidosus", 
"strPOTOs": "Poterioochromonas sp. DS", "strPSEUm": "Pseudo-nitzschia multistriata", 
"strPSEma": "Pseudo-nitzschia multistriata", "strPSEms": "Pseudo-nitzschia multiseries",
"strPTERd": "Pteridomonas danica", "strPYTHu": "Pythium ultimum var. sporangiiferum BR650", 
"strPYTHv": "Pythium vexans", "strRHIZC": "Rhizochromulina marina", "strSAPRd": "Saprolegnia diclina", 
"strSAPRp": "Saprolegnia parasitica", "strSCHYa": "Schizochytrium aggregatum",
"strSKEco": "Skeletonema costatum", "strSPUMv": "Spumella vulgaris", 
"strSTAUR": "Staurosira sp.", "strSYNCR": "Synchroma pusillum", "strTHAoc": "Thalassiosira oceanica",
"strTHAps": "Thalassiosira pseudonana", "strTHNEM": "Thalassionema frauenfeldii", 
"strTHNEn": "Thalassionema nitzschioides", "strTHRAc": "Thraustotheca clavata", "strTHRAU": "Thraustochytrium sp.", 
"strTHTRX": "Thalassiothrix antarctica", "strVAUCl": "Vaucheria litorea"}


basecolours = {'blue': '0000ff', 'brown': '996633', 'cyan': '00ffff', 'green': '00ff00', 'magenta': 'ff00ff', 'orange': 'ff8000', 'purple': '800080', 'red': 'ff0000' , 'yellow': 'ffff00', 'white': 'ffffff'}
black = ['-16777216', '000000']
if filtercolour in basecolours:
	filtercolour = basecolours[filtercolour]
elif filtercolour == 'all':
	print("any colour accepted")
else:
	print("unknown filter, setting to 'user-defined'. taxa with unrecognized colour codes will be retained")

#load fasta
seq_d = {}
badchars = ("|#=+,:;()'[]/") #also []/
for sequence in indataset:
	shortname = sequence.description.replace(" ","_")
	newname = []
	for c in shortname:
		if c in badchars:
			c = "_"
		newname.append(c)
	shortname = ''.join(newname)
	#print(shortname)
	seq_d[shortname] = sequence.seq
	if shortname.split("_")[0] in taxarepl9:
		taxon = taxarepl9[shortname.split("_")[0]]
		shortname = shortname.replace(shortname.split("_")[0], taxon.replace(" ", "_"))
		seq_d[shortname] = sequence.seq
	#print(">" + shortname + "\n" + seq_d[shortname])

print("done loading sequences")
#load taxa from tree
alltaxa = [] # a list of leaf labels with their color
pattern = r"\[&!color=.+"
skip = []
skippedc = 0
keptc = 0
treelines = intree.split('\n')[4:] #extract only leaf labels
for line in treelines:
	if line != ';':
		line = line.replace("'", "")
		line = line.replace("\t", "").replace(" ","_").replace("|","_")#.replace("@","_")
		if "%aligned" in line:
			line = clear_tag(line)
		#line = line.split("@")[0]
		#linecolour = line + colour
		#print(linecolour)
		alltaxa.append(line)
		if "[&!color" in line:
			colour = re.search(pattern, line).group()
			newcolour = line.split('[&!color=#')[1].replace("]", "")
			#print(taxon)
			if newcolour in black:
				print("black detected for %s, keeping this taxon" % (line))
				keptc += 1
			elif filtercolour == 'all':
				skip.append(line.split('[&!color=#')[0])
				skippedc += 1
			elif newcolour == filtercolour:
				skip.append(line.split('[&!color=#')[0])
				skippedc += 1
			else:
				print("unknown colour detected for %s, keeping this taxon" % (line))
				keptc += 1
			newtaxon = line.split('[&!color=#')[0]
			#print(newtaxon)
			alltaxa = [newtaxon if x==line else x for x in alltaxa] 
		else:
			#print(line)
			keptc += 1
			colour = ""

	else:
		break

errorfile = open("_key_errors.txt", "a")
print("done loading taxa")
if not args.no_omitted:
	print("omitted taxa listed in omitted-{}".format(args.fastain))
	#write omitted taxa
	with open('omitted-' + args.fastain, 'w') as f:
		for taxon in skip:
			nohighertaxon = taxon.split("@")[0]
			nospacetaxon = taxon.replace(" ","_")
			if seq_d.get(taxon) != None:
				f.write(">%s\n%s\n" % (taxon, seq_d[taxon]))
			elif seq_d.get(nohighertaxon) != None:
				f.write(">%s\n%s\n" % (taxon, seq_d[nohighertaxon]))
			elif seq_d.get(nospacetaxon) != None:
				f.write(">%s\n%s\n" % (taxon, seq_d[nospacetaxon]))
			else:
				print("!!!!!KEY ERROR for omitted", taxon)
				errorfile.write("omitted:\t{}\n".format(taxon))

print("writing filtered dataset...")
#write results
with open('{}-{}'.format(args.prefix, args.fastain), 'w') as out:
	for taxon in alltaxa:
		nohighertaxon = taxon.split("@")[0]
		nospacetaxon = taxon.replace(" ","_")
		if taxon not in skip:
			if seq_d.get(taxon) != None:
				out.write(">%s\n%s\n" % (taxon, seq_d[taxon]))
			elif seq_d.get(nohighertaxon) != None:
				out.write(">%s\n%s\n" % (taxon, seq_d[nohighertaxon]))
			elif seq_d.get(nospacetaxon) != None:
				out.write(">%s\n%s\n" % (taxon, seq_d[nospacetaxon]))
			else:
				print("!!!!!KEY ERROR for filtered", taxon)
				errorfile.write("filtered:\t{}\n".format(taxon))

print("WRITING DONE, \n\t{} taxa kept in {}-{},".format(keptc, args.prefix, args.fastain))
if not args.no_omitted:
	print("\t{} taxa omitted in omitted-{}".format(skippedc, args.fastain))
else:
	print("\t{} taxa omitted".format(skippedc))
print("Filtering finished!")
print("\n--------------------------------------------------------------------------")
