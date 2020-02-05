#!/usr/bin/python
#-*-coding:utf-8-*-

import random
from Bio.Seq import Seq
from statistics import *
import numpy as np
import sys
#seq=Seq(i)



for i in sys.argv:
	if i =="-R":
		R=int(sys.argv[sys.argv.index(i)+1])
	if i =="-out":
		out=sys.argv[sys.argv.index(i)+1]

tempo_out=out+"tempo/"
result=tempo_out+"data.fasta" #read split
source=tempo_out+"read_taille.fasta" #multifasta de read source

result=open(result,"w")
reads=open(source,"r")
tout_read=reads.readlines()

if len(tout_read) > 1000:
	R=500
elif len(tout_read) > 600 and len(tout_read) <1000:
	R=300
else:
	R=int(len(tout_read)/2)

liste_read=list()
for i in range(0,R):
	liste_read.append(random.randrange(0, len(tout_read)-2,2))
for l in liste_read:
	result.write(tout_read[l])
	result.write(tout_read[l+1])
print(len(tout_read),R)


result.close()
reads.close()
