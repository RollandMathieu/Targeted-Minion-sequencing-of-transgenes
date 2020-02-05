#!/usr/bin/python
#-*-coding:utf-8-*-

import random
from Bio.Seq import Seq
from statistics import *
import numpy as np
import sys
#seq=Seq(i)


Q=0
T=50
for i in sys.argv:
	if i =="-Q":
		Q=int(sys.argv[sys.argv.index(i)+1])
	if i =="-T":
		T=int(sys.argv[sys.argv.index(i)+1])
	if i =="-out":
		out=sys.argv[sys.argv.index(i)+1]

tempo_out=out+"tempo/"
result=tempo_out+"read_taille.fasta" #read split
source=tempo_out+"readsurcible.fasta" #multifasta de read source



result=open(result,"w")
reads=open(source,"r")
tout_read=reads.readlines()


seq=list()
for i in tout_read:
	if i.startswith(">"):
		seq.append(i[:-1])

Qscore=list()
taille=list()

for i in seq:
	ligne=i.split("_")
	Qscore.append(int(ligne[1]))
	taille.append(int(ligne[2]))

taillemin=np.percentile(taille,T)
Qscoremin=np.percentile(Qscore,Q)
#print(Qscoremin,taillemin)

compteur=0
while compteur < len(tout_read):
	if tout_read[compteur].startswith(">"):
		tempo=tout_read[compteur].split("_")
		if int(tempo[2]) > taillemin :#and int(tempo[1])>Qscoremin:
			result.write(tout_read[compteur]+tout_read[compteur+1])
	compteur=compteur+1

result.close()
reads.close()

#result=open(result2,"w")
#reads=open(result,"r")
#tout_read=reads.readlines()
