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

chrom = sys.argv[1]

featurelist = []
inversionlist = []

dembasfile = "DembeckBastideHits_positionsAlleles.txt"

inversionfile = "sixMajorInversions_rel5Coords.txt"

with open(inversionfile) as file:
	#ID,NAME,CHROM,BK1,BK2
	count = 0
	for line in file:
		if count>0:
			l = line.strip("\n").split("\t")
			if l[2]==chrom:
				entry = [l[0],l[1],l[2],int(l[3]),int(l[4])]
				inversionlist.append(entry)
		count+=1

gwashits = []
with open(dembasfile,"r") as file:
	#CHROM,POS,REF,ALT,SNP,DATASET
	count = 0
	for line in file:
		if count>0:
			l = line.strip("\n").split("\t")
			if l[0]==chrom:
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
			if (elements[FEATURE] in featureType) and (elements[CHROM]==chrom):
				entry = [elements[FEATURE],elements[CHROM],int(elements[START]),int(elements[END]),elements[STRAND]]
				featurelist.append(entry)

		except:
			continue

#same logic as other matching algorithm script, "DESTSNPs_MatchingAlgorithm_RegionInversionFreq_bootstrapSamps.py"
orchfile = "Orch16_ThreeTimepoint_Freqs_"+chrom+".txt"

orchdict = {}
regionPriority = {"exon":1,"three_prime_UTR":2,"five_prime_UTR":2,"intron":3,"gene":4,"intergenic":5}

with open(orchfile) as file:
	#ROWNUM,CHROM,POS,TP1,TP2,TP3
	count = 0
	for line in file:
		if count>0:
			l = line.strip("\n").split(" ")
			key = (l[1],int(l[2]))
			pos = int(l[2])

			region = "intergenic"
			for f in featurelist:
				if f[1]==chrom:
					start = f[2]
					end = f[3]

					if pos>=start and pos<=end:
						if f[0] != region:
							if regionPriority[f[0]]<regionPriority[region]:
								region = f[0]
			if chrom!="X":
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
			else:
				distance = "N/A"

			avgfreq = (float(l[3])+float(l[4])+float(l[5]))/3

			entry = [region,distance,avgfreq]
			orchdict.update({key:entry})
		count+=1

buffer = 50000
wiggle = 0.01

def matchedgroup(snp,snpdict,wiggle,buffer):
	matchedlist = []
	chrom = snp[0]
	pos = snp[1]

	hitinfo = snpdict[snp]
	region = hitinfo[0]
	distance = hitinfo[1]
	avgfreq = hitinfo[2]

	if distance=="N/A":
		invertstatus = "N/A"
	elif distance==0:
		invertstatus = "inside"
	elif distance<500000:
		invertstatus = "close"
	else:
		invertstatus = "far"

	highfreq = avgfreq+wiggle
	lowfreq = avgfreq-wiggle

	if lowfreq<0:
		lowfreq=0


	for s in snpdict.keys():
		if (s[0]==chrom) and (abs(s[1]-pos)>=buffer):
			sinfo = snpdict[s]
			if sinfo[0]==region:
				d = sinfo[1]
				if d=="N/A":
					i = "N/A"
				elif d==0:
					i = "inside"
				elif d<500000:
					i = "close"
				else:
					i = "far"

				if i==invertstatus:
					if (sinfo[2]>=lowfreq) and (sinfo[2]<=highfreq):
						matchedlist.append(s)

	if len(matchedlist)>=200:
		return(matchedlist)
	else:
		wiggle = wiggle+0.01
		return(matchedgroup(snp,snpdict,wiggle,buffer))



matchsnpfile = "MatchedSNPs_100Groups_Orch16_"+chrom+".txt"
groupsnpfile = "GroupPermutations_Orch16_"+chrom+".txt"


with open(matchsnpfile,"w") as matchfile:
	matchfile.write("SAMP\tHITCHROM\tHITPOS\tMATCHCHROM\tMATCHPOS\n")
	with open(groupsnpfile,"w") as groupfile:
		groupfile.write("GROUP\tHITCHROM\tHITPOS\tMATCHCHROM\tMATCHPOS\n")
		for hit in gwashits:
			snp = (hit[0],hit[1])
			if snp in orchdict.keys():
				matches = matchedgroup(snp,orchdict,wiggle,buffer)

				for i in range(1,101):
					samp = random.sample(matches,100)
					for s in samp:
						matchfile.write(str(i)+"\t"+str(snp[0])+"\t"+str(snp[1])+"\t"+str(s[0])+"\t"+str(s[1])+"\n")

				for i in range(1,1001):
					perm = random.choice(matches)
					groupfile.write(str(i)+"\t"+str(snp[0])+"\t"+str(snp[1])+"\t"+str(perm[0])+"\t"+str(perm[1])+"\n")
			else:
				print("no match "+str(snp))



