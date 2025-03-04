from ftplib import FTP
import gzip
import re
import sys
import copy
import os
import subprocess
import math
import statistics
import random


## This code takes a series of candidate SNPs from two previous studies (Bastide et al. 2013 and Dembeck et al. 2015)
## and determines 1. the genetic region in which they're located (exon, 5'/3' UTR, intron, generally a "gene", or intergenic, 
## in order of priority), 2. their inversion status (inside one of six common inversions found in natural populations -- aka 
## a distance of "0", or the distance to the nearest inversion breakpoint), and 3. the average frequency of the SNP across the 5 DEST samples that
## we're using for our clinal analysis. 

## The DEST file has to be subset to the clinal samples and converted into a text file before this code can be run -- DESTgenes_fivesamps.py

chroms = ["2L","2R","3L","3R","X"]



## This is a series of functions that convert release 5 dmel coordinates to release 6 dmel coordinates and vice versa.
## The file "Schroeder.Dmel_r5-to-r6_map.2014.6.12.txt" is from Flybase

mapfile = "SupplementaryFiles/Schroeder.Dmel_r5-to-r6_map.2014.6.12.txt"
mapswap = []

#reads in map file -- each entry is a region that has been shifted as a block to new release 6 coordinates
#each entry contains, starting at index 0 - rel 5 chromosomal arm, 1 - rel 5 BP start, 2 - rel 5 BP end,
#	3 - rel 6 chromosomal arm, 4 - rel 6 start, 5 - rel 6 end, 6 - rel 6 strand
with open(mapfile) as file:
	for line in file:
		if line.startswith("#"):
			continue
		else:
			l = line.strip("\n").split("\t")
			if l[0] in chroms:
				mapswap.append([l[0],int(l[1]),int(l[2]),l[3],int(l[4]),int(l[5]),l[6]])


#each of these functions checks in which region a SNP is located in either rel 5 or rel 6 coordinates and then
#shifts the given position by the mapping
def R5toR6(poschrom,maplist):
	pos = poschrom[1]
	chrom = poschrom[0]
	newchrom = ""

	for i in maplist:
		if i[0]==chrom:
			if pos>=i[1] and pos<=i[2]:
				newchrom = i[3]
				start = i[1]
				end = i[2]
				newstart = i[4]
				newend = i[5]
				if start==newstart and end==newend:
					newpos = pos
				else:
					dif = pos-start
					newpos = dif+newstart
					
	if not newchrom:
		return "NA"
	else:
		return (newchrom,newpos)

def R6toR5(poschrom,maplist):
	pos = poschrom[1]
	chrom = poschrom[0]
	newchrom = ""

	for i in maplist:
		if i[3]==chrom:
			if pos>=i[4] and pos<=i[5]:
				newchrom = i[0]
				start = i[4]
				end = i[5]
				newstart = i[1]
				newend = i[2]
				if start==newstart and end==newend:
					newpos = pos
				else:
					dif = pos-start
					newpos = dif+newstart
	if not newchrom:
		return "NA"
	else:
		return (newchrom,newpos)





#info for downloading gene annotations from Flybase. We're using a gene annotation file from rel 5 to keep things
#consistent for the Bastide/Dembeck hits, which are in release 5
ftpAddress = "ftp.flybase.org"
gtfDirectory = "releases/FB2014_03/dmel_r5.57/gff"
gtfFile = "dmel-all-r5.57.gff.gz"
featureType = ["exon","three_prime_UTR","five_prime_UTR","gene","intron"]

#gtf file positions
CHROM = 0
FEATURE = 2
START = 3
END = 4
STRAND = 6
GENEINFO = 8


#empty data structures
featurelist = []
inversionlist = []
gwashits = []

#a list of every SNP that reached significance in Bastide et al 2013 and Dembeck et al 2015 (in rel 5)
dembasfile = "SupplementaryFiles/DembeckBastideHits_positionsAlleles.txt"

#a list of the 6 common inversions found in natural populations, coordinates in rel 5
inversionfile = "SupplementaryFiles/sixMajorInversions_rel5Coords.txt"

#reads in inversion coordinates -- 0 - flybase ID, 1 - name, 2 - chromosomal arm, 3 - breakpoint 1 position, 4 - breakpoint 2 position
with open(inversionfile) as file:
	count = 0
	for line in file:
		if count>0:
			l = line.strip("\n").split("\t")
			print(l)
			entry = [l[0],l[1],l[2],int(l[3]),int(l[4])]
			inversionlist.append(entry)
		count+=1

#reads in Bastide/Dembeck GWAS hits -- 0 - chromosomal arm, 1 - rel 5 position, 2 - reference allele, 3 - alternative allele, 4 - light-associated
#	allele in Dembeck et al 2015, 5 - dataset that hit belongs to (either "Bastide" or "Dembeck")
with open(dembasfile,"r") as file:
	count = 0
	for line in file:
		if count>0:
			l = line.strip("\n").split("\t")
			entry = [l[0],int(l[1]),l[2],l[3],l[4],l[5]]
			gwashits.append(entry)
		count+=1



#download and read in flybase annotation file and get all genes, exon/intron/5' & 3' UTR positions
if not os.path.isdir("SupplementaryFiles/tempDmelFiles"):
	subprocess.call(["mkdir","tempDmelFiles"],shell=True)
	print("folder tempDmelFiles created")

#if it doesn't already exist, copy gtf file from flybank
if not os.path.isfile("SupplementaryFiles/tempDmelFiles/dmel2014.gff.gz"):
	ftp = FTP(ftpAddress)
	ftp.login()
	ftp.cwd(gtfDirectory)
	#make the right call string from the gtfFile created
	retrstring = "RETR "+gtfFile
	ftp.retrbinary(retrstring,open("SupplementaryFiles/tempDmelFiles/dmel2014.gff.gz", 'wb').write)
	print("temp file dmel2014.gff.gz created from flybase gff file")

count = 0
#read in the gtf file
# 0 - feature type, 1 - chromosomal arm, 2 - feature start, 3 - feature end, 4 - feature strand
with gzip.open("SupplementaryFiles/tempDmelFiles/dmel2014.gff.gz","rt") as file:
	for line in file:
		if line.startswith("##"):
			continue
		elements = line.strip("\n").split("\t")
		#select for elements that are the right feature type
		try:
			if elements[FEATURE] in featureType:
				entry = [elements[FEATURE],elements[CHROM],int(elements[START]),int(elements[END]),elements[STRAND]]
				featurelist.append(entry)

		except:
			continue

freqdict = {}
unsorted = []
nomap = []

#open DEST frequencies for the 5 clinal samples and convert to rel 5 (currently in rel 6)
#calculate the mean frequency across the 5 clinal DEST samples & record the number of samples that the SNP appears in
for chrom in chroms:
	chromdata = []
	filename = "DESTdata_"+chrom+"_pooledsampsfreq_fivesamps.txt"
	with open(filename, "r") as file:
		count = 0
		for line in file:
			if count>0:
				l = line.strip("\n").split("\t")
				chrompos = [l[0],int(l[1])]
				rel5chrompos = R6toR5(chrompos,mapswap)
				if rel5chrompos !="NA":
					actualchrom = rel5chrompos[0]
					actualpos = int(rel5chrompos[1])
					if actualchrom != chrom:
						unsorted.append(l)
						continue
					entry = [actualchrom,actualpos,l[2],l[3]]
					existingfreq = []
					numsamp = 0
					for i in range(4,len(l)):
						if l[i] != ".":
							entry.append(float(l[i]))
							existingfreq.append(float(l[i]))
							numsamp+=1
						else:
							entry.append(l[i])
					if existingfreq:
						entry.append(statistics.mean(existingfreq))
					else:
						entry.append("N/A")
					entry.append(numsamp)
					chromdata.append(entry)
				else:
					nomap.append(l)

			count+=1

	freqdict.update({chrom:chromdata})

print(len(unsorted))
print(len(chromdata))
print(len(nomap))







regionPriority = {"exon":1,"three_prime_UTR":2,"five_prime_UTR":2,"intron":3,"gene":4,"intergenic":5}
gwasprint = []

#determine which type of genetic region the SNP is located in -- the highest "priority" label is applied to the SNP
#so if a SNP is in the exon of one gene and the intron of another, it gets the "exon" label
for snp in gwashits:
	chrom = snp[0]
	pos = snp[1]

	ref = snp[2]
	alt = snp[3]
	snpass = snp[4]
	dataset = snp[5]

	freqdata = freqdict[chrom]

	regions = []
	inversiondists = []
	region = "intergenic"

	for f in featurelist:
		if f[1]==chrom:
			start = f[2]
			end = f[3]

			if pos>=start and pos<=end:
				regions.append(f)
				if f[0] != region:
					if regionPriority[f[0]]<regionPriority[region]:
						region = f[0]

#determine the distance to the nearest common inversion
	for i in inversionlist:
		if i[2]==chrom:
			start = i[3]
			end = i[4]
			distance = 9000000000

			if pos>=start and pos<=end:
				distance = 0
			else:
				startdis = abs(start-pos)
				enddis = abs(end-pos)
				if startdis<distance:
					distance = startdis
				if enddis<distance:
					distance = enddis
			inversiondists.append(distance)

	numallele=0
	inDEST=False
	inAlleleDEST = False

	for f in freqdata:
		if f[0]==chrom and f[1]==pos:
			inDEST = True
			if f[2]!=ref:
				print("ref issue")
				print(snp,f)
			else:
				if f[2]==ref and f[3]==alt:
					inAlleleDEST=True
					meanfreq = f[-2]
					numDESTsamps = f[-1]

	if inDEST and not inAlleleDEST:
		print("allele mismatch")
		print(snp)

	if inAlleleDEST:
		entry = [chrom,pos,ref,alt,snpass,dataset,region,meanfreq,numDESTsamps]
	else:
		entry = [chrom,pos,ref,alt,snpass,dataset,region,"N/A",0]

	if inversiondists:
		entry.append(min(inversiondists))
	else:
		entry.append("N/A")

	gwasprint.append(entry)




outputfile = "GWAShits_DembeckBastide_ClinalDEST_ForMatching.txt"

with open(outputfile, "w") as wfile:
	wfile.write("CHROM\tPOS\tREF\tALT\tSNPASSOC\tDATASET\tREGION\tFREQ\tNUMDESTSAMPS\tINVERSIONDIST\n")
	for i in gwasprint:
		wfile.write(str(i[0]))
		for j in range(1, len(i)):
			wfile.write("\t"+str(i[j]))
		wfile.write("\n")