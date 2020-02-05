#!/usr/bin/python
#-*-coding:utf-8-*-

import sys
from Bio.Seq import Seq

for i in sys.argv:
	if i =="-out":
		out=sys.argv[sys.argv.index(i)+1]


tempo_out=out+"tempo/"
t=tempo_out+"read.fasta" # read
result=tempo_out+"readsurcible.fasta" #read sélectionner
source=tempo_out+"multifasta.fasta" #multifasta de source



result=open(result,"w")
reads=open(t,"r")
sources=open(source,"r")

tout_read=reads.readlines()
tout_source=sources.readlines()

#tout_read=list(set(tout_read2))

nb=0 #nombre de read qui match
nj=len(tout_source)/2 #nombre de read total

for l in tout_read:
	ligne2=l.split(",")
	numread=ligne2[1]
	i=str(numread)
	y=0
	while y <len(tout_source):
		ligne=tout_source[y]
		if i[:-1]==ligne[1:-1]:
			nb=nb+1
			if ligne2[0] =="plus":
				result.write(">"+i+tout_source[y+1])
				y=len(tout_source)
			else :
				seq=Seq(tout_source[y+1])
				result.write(">"+i[:-1]+str(seq.reverse_complement())+"\n")
				y=len(tout_source)
		y=y+1

print(nb," reads qui match sur la cible sur",nj,"soit ",nb/(nj)*100,"%")

result.close()
reads.close()
sources.close()
