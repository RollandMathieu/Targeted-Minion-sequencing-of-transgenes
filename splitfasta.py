#!/usr/bin/python
#-*-coding:utf-8-*-

import random
from Bio.Seq import Seq
import sys

#seq=Seq(i)

for i in sys.argv:
	if i =="-out":
		out=sys.argv[sys.argv.index(i)+1]

tempo_out=out+"tempo/"
result=tempo_out+"readsplit.fasta" #read split
source=tempo_out+"data.fasta" #multifasta de read source
#source="readsurcible.fasta"

result=open(result,"w")
reads=open(source,"r")

tout_read=reads.readlines()

nb=1 #nombre de read
R=200 #taille des read voulu

for i in tout_read:
	if i.startswith(">") is False:
		DNA=i[:-1]
		#nr=len(dna)//R
		p=0
		while p<len(DNA)-R:
			result.write(">"+str(nb)+"\n"+DNA[p:p+R]+"\n")
			nb=nb+1
			p=p+R
		result.write(">"+str(nb)+"\n"+DNA[p:]+"\n")
		nb=nb+1
result.close()
reads.close()
