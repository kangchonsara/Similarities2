#Calculate amino acid similarity of sequences at each year to 2016-2017 vaccine sequence at site A and site B
import os

beginY = 1968

def read_epitope_shih(epifName):
	sites = ['a','b','c','d']
	epitopes = [[] for i in range(len(sites)) ]
	
	epif = open(epifName, "r")
	for line in epif:
		each = line.split("\n")[0].split("\t")
		idx = sites.index(each[1])
		epitopes[idx].append(int(each[0]))
			
	epif.close()
	
	return epitopes

def read_epitope_koel(epifName):
	
	epitope = []
	
	epif = open(epifName, "rU")
	for line in epif:
		epitope.append(int(line.split("\n")[0]))
		
	epif.close()
	
	return epitope
	
def read_vaccineSeq(vacfName):
	vacSeqs = []
	vacf = open(vacfName, "r")
	for line in vacf:
		if line.find(">") >= 0:
			seq = line.split("\n")[0] + " "
		else:
			seq += line.split("\n")[0]
			vacSeqs.append(seq)
			
	vacf.close()
	return vacSeqs
	
def vacSeq_dif(vacSeqs, epitope):
	vs_dif_all = 0
	num_comp = 0
	for v1 in range(len(vacSeqs)-1):
		for v2 in range(v1+1, len(vacSeqs)):
			vac1 = vacSeqs[v1].split(" ")[1]
			vac2 = vacSeqs[v2].split(" ")[1]
			
			num_comp += 1
			vs_dif = 0
			for s in epitope:
				if vac1[s-1] != vac2[s-1]:
					print (s, v1, v2, vac1[s-1], vac2[s-2])
					vs_dif += 1

			vs_dif_all += vs_dif
	vs_dif_avg = 1.0*vs_dif_all/(len(epitope)*num_comp)
	
	return vs_dif_avg
	
			
	
#When there are 2 or more identical samples (same virus name and same sequence), leave only 1.
#Remove samples with deletion or ambiguity at epitope site.
def remove_dupl_ambig(sequences, epitope_mixed):
	years = []

	for year in sequences:

		oneYear = []
		for seq in year:
			if seq in oneYear:
				continue
			else:
				sequence = seq.split(" ")[1]
				
				remove = 0
				for s in range(len(sequence)):
					if s+1 in epitope_mixed:
						if sequence[s+1] == "-" or sequence[s+1] == "?" or sequence[s+1] == "*":
							remove = 1
							break
							
				if remove == 1:
					continue

				oneYear.append(seq)
				
		years.append(oneYear)
	return years

def get_seqs_byYear(infName, beginY, endY):
	inf = open(infName, "r")

	AllbyYear_r = [[] for i in range(beginY, endY+1)]
	USbyYear_r = [[] for i in range(beginY, endY+1)]
	for line in inf:
		if line.find(">") >= 0:
			year = int(line.split("|")[2].split("/")[0])
			country = line.split("|")[3]
			seq = line.split("|")[1] + " "
		else:
			if year < beginY or year > endY:
				continue
			seq += line.split("\n")[0]
			AllbyYear_r[year-beginY].append(seq)
			if country == "USA":
				USbyYear_r[year-beginY].append(seq)
	inf.close()
	
	return AllbyYear_r, USbyYear_r

def calc_epi_similarity(seqs_byYear, vacSeqs, epitope):
	epi_similarities = []
	
	for yearSeq in seqs_byYear:
		mean_simil = 0
		
		for seq in yearSeq:
			seq = seq.split(" ")[1]
			similarity = 0
			for vacSeq in vacSeqs:
				vacSeq = vacSeq.split(" ")[1]
				for s in epitope:
					if seq[s-1] == vacSeq[s-1]:
						similarity += 1
			
			similarity = 1.0*similarity/(len(epitope)*len(vacSeqs)) #per site
			mean_simil += similarity
			
		mean_simil = mean_simil/len(yearSeq) #per seq
		epi_similarities.append(mean_simil)
	
	return epi_similarities

def calc_gly_proportion(seqs_byYear, vacSeqs, gly_site, firstY):

	s = gly_site[0]
	allele = gly_site[1]
	startY_idx = gly_site[2] - firstY
	endY_idx = gly_site[3] - firstY
	
	gly_proportions = []
	for y in range(len(seqs_byYear)):
		proportion = 0

		if y < startY_idx or y > endY_idx:
			gly_proportions.append(proportion)
			continue
		
		for seq in seqs_byYear[y]:
			seq = seq.split(" ")[1]
			if seq[s-1] == allele:
				proportion += 1
				
		proportion = 1.0*proportion/len(seqs_byYear[y])
		gly_proportions.append(proportion)
		
	return gly_proportions
			
			
	
	
#parameters
usingKseqs = 'Any' #'Any' or 'K'
vac = 'Egg' #'Egg', 'HK', "SW"
firstY = 1968 #first year of analysis
lastY = 2016 #last year of analysis
exclusions = [160]

#input files
shifName = os.path.normpath("../data/shih_epitope.txt")
koel_epitope_bfName = os.path.normpath("../data/koel_epitope_b.txt")
koel_cluster_changing_bfName = os.path.normpath("../data/koel_cluster_changing_b.txt")

if vac == "Egg":
	vacfName = os.path.normpath("../data/egg_adapted_vaccine_20162017_AA.fasta")
elif vac == "HK":
	vacfName = os.path.normpath("../data/HongKong4801_AA.fas")
elif vac == "SW":
	vacfName = os.path.normpath()
	
#output file
similfName = os.path.normpath("../result/similarities_")
glyfName = os.path.normpath("../result/proportion_gly_")
AllbyYear_r = [] 
USbyYear_r =[]

#from 1968 to 2010, read from ncbi
beginY = 1968
endY = 2010
infName = os.path.normpath("../data/H3_6810AA.fasta")
AllbyYear_temp, USbyYear_temp = get_seqs_byYear(infName, beginY, endY)
AllbyYear_r += AllbyYear_temp
USbyYear_r += USbyYear_temp

#from 2011 to 2011, read from ncbi
beginY = 2011
endY = 2011
infName = os.path.normpath("../data/ncbi_aligned_6812_AA.fas")
AllbyYear_temp, USbyYear_temp = get_seqs_byYear(infName, beginY, endY)
AllbyYear_r += AllbyYear_temp
USbyYear_r += USbyYear_temp

#from 2012 to 2016, read from gisaid
beginY = 2012
endY = 2016
infName = os.path.normpath("../data/gisaid_aligned_1217_AA.fas")
AllbyYear_temp, USbyYear_temp = get_seqs_byYear(infName, beginY, endY)
AllbyYear_r += AllbyYear_temp
USbyYear_r += USbyYear_temp

#read epitope sites 
epitopes_shih = read_epitope_shih(shifName)
koel_epitope_b = read_epitope_koel(koel_epitope_bfName)
koel_cluster_changing_b = read_epitope_koel(koel_cluster_changing_bfName)
epitope_variants = [epitopes_shih[0], epitopes_shih[1], koel_epitope_b, koel_cluster_changing_b]
epitope_names = ["siteA_shih", "siteB_shih", "siteB_koel_clusterDifference", "siteB_koel_clusterChanging"]
epitope_mixed = []
for epitope in epitope_variants:
	epitope_mixed += epitope
epitope_mixed = list(set(epitope_mixed))

#remove sites in exclusions from epitopes
for e in range(len(epitope_variants[:])):
	for exc in exclusions:
		try:
			epitope_variants[e].remove(exc)
		except ValueError:
			continue


#glycosylated site
gly_sites = [(160, 'T', 1978, 2017)]

#remove duplicated sequences or sequences with ambiguous allele in epitope sites  
AllbyYear = remove_dupl_ambig(AllbyYear_r, epitope_mixed)
USbyYear = remove_dupl_ambig(USbyYear_r, epitope_mixed)

#for year with number of US sequences >= 70, us US sequences. 
#otherwise, use all sequences available
isUSseq = []
seqs_byYear = []
for y in range(len(USbyYear)):
	if len(USbyYear[y]) >= 70:
		year = USbyYear[y]
		isUSseq.append(1)
	else:
		year = AllbyYear[y]
		isUSseq.append(0)
	seqs_byYear.append(year)

	
#When using sequences having K in 160 only
Kseqs_byYear = []
for y in range(1968, endY+1):
	Kseqs_1year = []
	ydx = y-beginY
	for seq in seqs_byYear[ydx]:
		if seq.split(" ")[1][160-1] == 'K':
			Kseqs_1year.append(seq)
	Kseqs_byYear.append(Kseqs_1year)
	
if usingKseqs == 'K':
	seqs_byYear = Kseqs_byYear

#read vaccine sequence	
vacSeqs = read_vaccineSeq(vacfName)
vacSeqs = vacSeqs[0:1]

#Now calculated similarity

for e in range(len(epitope_variants)):
	similf = open(similfName+epitope_names[e]+".csv", "w")
	similf.write("year,similarity,numSites\n")
	similarities = calc_epi_similarity(seqs_byYear, vacSeqs, epitope_variants[e])
	for y in range(len(similarities)):
		if usingKseqs == 'K':
			firstY = 1978
		similf.write(str(y+firstY)+","+str(similarities[y])+","+str(len(epitope_variants[e]))+"\n")
	similf.close()

for e in range(len(exclusions)):
	similf = open(similfName+str(exclusions[e])+".csv", "w")
	similf.write("year,similarity\n")
	similarities = calc_epi_similarity(seqs_byYear, vacSeqs, [exclusions[e]])
	for y in range(len(similarities)):
		similf.write(str(y+firstY) + "," + str(similarities[y])+"\n")
	similf.close()
	
for g in range(len(gly_sites)):
	glyf = open(glyfName+str(gly_sites[g][0])+".csv", "w")
	glyf.write("year,gly_proportion\n")
	proportions = calc_gly_proportion(seqs_byYear, vacSeqs, gly_sites[g], firstY)
	for y in range(len(proportions)):
		glyf.write(str(y+firstY) + "," + str(proportions[y])+"\n")
	glyf.close()
		
	













	
	
	
	
	
	
	