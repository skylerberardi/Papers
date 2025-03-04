#code to pull existing hits (from Bastide et al. 2013 and Dembeck et al. 2015) from Linvilla, PA seasonal DEST samples
#for each Dembeck/Bastide hit, convert the dmel genome rel 5 coordinates to rel 6, determine the inversion status, the genetic region (exon, intron, etc.), and the avg frequency in the seasonal samples 

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
import pysam
from scipy.stats import hypergeom
from scipy.stats import kstest
from scipy.stats import ks_2samp
from sympy import symbols
from sympy.solvers import solve

#these are the linvilla samples
pooledsampsID = ["PA_li_09_spring","PA_li_09_fall","PA_li_10_spring","PA_li_10_fall","PA_li_11_spring","PA_li_11_fall","PA_li_12_spring","PA_li_12_fall","PA_li_14_spring","PA_li_14_fall","PA_li_15_spring","PA_li_15_fall"]
pooledsampsdict = {"PA_li_09_spring":[9,"spring",2009.0],"PA_li_09_fall":[9,"fall",2009.5],"PA_li_10_spring":[10,"spring",2010.0],"PA_li_10_fall":[10,"fall",2010.5],"PA_li_11_spring":[11,"spring",2011.0],"PA_li_11_fall":[11,"fall",2011.5],"PA_li_12_spring":[12,"spring",2012.0],"PA_li_12_fall":[12,"fall",2012.5],"PA_li_14_spring":[14,"spring",2014.0],"PA_li_14_fall":[14,"fall",2014.5],"PA_li_15_spring":[15,"spring",2015.0],"PA_li_15_fall":[15,"fall",2015.5]}

chroms = ["2L","2R","3L","3R","X"]


mapfile = "Schroeder.Dmel_r5-to-r6_map.2014.6.12.txt"
mapswap = []

with open(mapfile) as file:
	for line in file:
		if line.startswith("#"):
			continue
		else:
			l = line.strip("\n").split("\t")
			if l[0] in chroms:
				mapswap.append([l[0],int(l[1]),int(l[2]),l[3],int(l[4]),int(l[5]),l[6]])

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



featurelist = []
inversionlist = []

dembasfile = "DembeckBastideHits_positionsAlleles.txt"

inversionfile = "sixMajorInversions_rel5Coords.txt"

with open(inversionfile) as file:
	count = 0
	for line in file:
		if count>0:
			l = line.strip("\n").split("\t")
			print(l)
			entry = [l[0],l[1],l[2],int(l[3]),int(l[4])]
			inversionlist.append(entry)
		count+=1

gwashits = []
with open(dembasfile,"r") as file:
	count = 0
	for line in file:
		if count>0:
			l = line.strip("\n").split("\t")
			entry = [l[0],int(l[1]),l[2],l[3],l[4],l[5]]
			gwashits.append(entry)
		count+=1

#if there doesn't exist a temp file for files created during this -
if not os.path.isdir("tempDmelFiles"):
	#then create it
	#this isn't a option in the config file to prevent code injection
	#that's very unlikely, but better safe than sorry
	subprocess.call(["mkdir","tempDmelFiles"],shell=True)
	print("folder tempDmelFiles created")

#if it doesn't already exist, copy gtf file from flybank
if not os.path.isfile("tempDmelFiles/dmel2014.gff.gz"):
	ftp = FTP(ftpAddress)
	ftp.login()
	ftp.cwd(gtfDirectory)
	#make the right call string from the gtfFile created
	retrstring = "RETR "+gtfFile
	ftp.retrbinary(retrstring,open("tempDmelFiles/dmel2014.gff.gz", 'wb').write)
	print("temp file dmel2014.gff.gz created from flybase gff file")

count = 0
#read in the gtf file
with gzip.open("tempDmelFiles/dmel2014.gff.gz","rt") as file:
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
for chrom in chroms:
	chromdata = []
	filename = "DESTdata_"+chrom+"_pooledsampsfreq_LinvillaSeason.txt"
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
						sampidindex = i-4
						if l[i] != ".":
							sampid = pooledsampsID[sampidindex]
							sampseason = pooledsampsdict[sampid][1]
							if sampseason=="fall":
								previousfreq = entry[-1]
								if previousfreq!=".":
									numsamp+=1
							entry.append(float(l[i]))
							existingfreq.append(float(l[i]))
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




outputfile = "GWAShits_DembeckBastide_SeasonDEST_ForMatching.txt"

with open(outputfile, "w") as wfile:
	wfile.write("CHROM\tPOS\tREF\tALT\tSNPASSOC\tDATASET\tREGION\tFREQ\tNUMDESTSAMPS\tINVERSIONDIST\n")
	for i in gwasprint:
		wfile.write(str(i[0]))
		for j in range(1, len(i)):
			wfile.write("\t"+str(i[j]))
		wfile.write("\n")