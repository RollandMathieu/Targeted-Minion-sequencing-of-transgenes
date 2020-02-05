#!/usr/bin/python
#-*-coding:utf-8-*-

import sys
from Bio.Seq import Seq

#seq=Seq(i)
for i in sys.argv:
	if i =="-out":
		out=sys.argv[sys.argv.index(i)+1]


tempo_out=out+"tempo/"

t=tempo_out+"scaffolds.fasta" # read
result=tempo_out+"scaffolds2.fasta" #read sélectionner

result=open(result,"w")
reads=open(t,"r")

tout_fichier2=reads.readlines()
tout_fichier=tout_fichier2[:-1]

compteur= 0

while compteur <len(tout_fichier):
	dna=str()
	ligne=tout_fichier[compteur]
	if ligne.startswith(">"):
		result.write(ligne)
		compteur=compteur+1
	else:
		while tout_fichier[compteur].startswith(">") is not True and compteur<len(tout_fichier)-1:
			ligne=tout_fichier[compteur]
			dna=dna+ligne[:-1]
			if compteur <= len(tout_fichier):
				compteur=compteur+1
		if compteur ==len(tout_fichier)-1:
			compteur=compteur+1
		result.write(dna+"\n")

result.close()
reads.close()
