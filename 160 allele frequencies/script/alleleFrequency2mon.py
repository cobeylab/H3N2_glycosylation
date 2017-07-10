#Calculate 160 allele frequencies from Dec 2011 to Mar 2017 in two-month-window

beginY = 2011
endY = 2017


def remove_dupl_ambig(sequences):
#remove same samples (same sequence and same name)
#and also remove samples with ambiguity or deletion at site 160
	years = []
	for year in sequences:
		oneYear = []
		for seq in year:
			if seq in oneYear:
				continue
			else:
				site = seq.split(" ")[1][159]
				if site == "-" or site == "?":
					continue
				oneYear.append(seq)
		years.append(oneYear)
	return years

def check_labstrains():
	fname = '_loc_age.txt'
	for y in range(2012, 2018):
		inf = open("../data/"+str(y)+fname, "r")
		for line in inf:
			if line.find("Location") >= 0: 
				each = line.split("\t")
				colLoc = each.index("Location")
				colHA = each.index("HA Segment_Id")
			else:
				each = line.split("\t")			
				if each[colLoc].lower().find("lab") == 1:
					print (each[colLoc])

		inf.close()

def byYear_gisaid(fastaf):
	inf = open(fastaf, "r")
	byYear = [[] for i in range(beginY, endY+1)]
	for line in inf:
		if line.find(">") >= 0:
			try:
				year = int(line.split("_|_")[2].split("-")[0])
			except ValueError:
				year = int(line.split("_|_")[2].split("_")[0])
			seq = line.split("_|_")[1] + "_|_" + line.split("_|_")[2] + " "
		else:
			seq += line.split("\n")[0]
			byYear[year-beginY].append(seq)

	inf.close()
	return byYear
	
def byMonth(year):
	months = [[] for i in range(12)]
	for seq in year:
		time = seq.split(" ")[0].split("_|_")[1]
		if len(time.split("-")) < 2:
			continue
		try:
			month =int(time.split("-")[1])
		except ValueError:
			month = int(time.split("-")[1].split("_")[0])
		months[month-1].append(seq)
	return months

def by2month(oneMon):
	twoMon = []
	for m in range(0, len(oneMon), 2):
		twoMon.append(oneMon[m]+oneMon[m+1])
	return twoMon
	
def calc_frequency(byYear):
	frequency = []
	for i in range(len(byYear)):
		if len(byYear[i]) == 0:
			continue
		alleles = [0,0,0]
			
		for virus in byYear[i]:
			seq = virus.split(" ")[1]
			if seq[159] == "T": 
				alleles[0] += 1
			elif seq[159] == "K":
				alleles[1] += 1
			else:
				alleles[2] += 1

		s = sum(alleles)
		freq = []
		for a in alleles:
			f = round( a/s, 3)
			freq.append(f)
		frequency.append(freq)
		
	return frequency

	
#confirm that there is no lab strains
check_labstrains()

#put sequences into a list by year first
fastaf = "../data/20112017_align_AA.fasta"
byYear_r = byYear_gisaid(fastaf)
byYear = remove_dupl_ambig(byYear_r)

#put sequences into a list by month
months = [] 
for y in range(len(byYear)):
	months += byMonth(byYear[y])

#put sequences into a list by two-month-window, 
#from Dec 2011 to Mar 2017
months = by2month(months[11:11+64]) #indexing for: from Dec 2011 to Mar 2017
frequency = calc_frequency(months)


for idx in range(len(frequency)):
	i = frequency[idx]
	print (str(i[0]) + "," + str(i[1]))

	
	
	
	
	
	
	
	
	
	
	
	
	
	