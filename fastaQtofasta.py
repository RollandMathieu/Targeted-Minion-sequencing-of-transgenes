#!/usr/bin/python
#-*-coding:utf-8-*-
import os
from Bio.Seq import Seq
from Bio import SeqIO
from statistics import median
import sys
# récupérer les fasta d'un fastq si plus long que n  et ajoute le reverse complément de chaque read

n=3000
for i in sys.argv:
	if i =="-reads":
		reads=sys.argv[sys.argv.index(i)+1]
	if i =="-out":
		out=sys.argv[sys.argv.index(i)+1]
	if i =="-len":
		n=int(sys.argv[sys.argv.index(i)+1])
tempo_out=out+"tempo/"
merge=os.listdir(reads)


nb=1 #nombre de reads sauvegarders
nj=0 #nombre de reads jeter
t=tempo_out+"fastq.fastq" # input
multifasta=tempo_out+"multifasta.fasta" #output
multifast=open(multifasta,"w")
multifast.close()


for y in merge:
	ex=y.split(".")
	if ex[-1] =="fastq":
		#fichier=open(y,"r") #source
		multifast=open(multifasta,"a")
		#tout_fichier=fichier.readlines()
		chemin=reads+y
		for rec in SeqIO.parse(chemin, "fastq"):
			if  len(rec.seq)>=n:
				multifast.write(">"+str(nb)+"_"+str(int(median(rec.letter_annotations["phred_quality"])))+"_"+str(len(rec.seq))+"\n"+str(rec.seq)+"\n")
				nb=nb+1
			else:
				nj=nj+1
		#fichier.close()
		multifast.close()
print(nb-1," reads sauvegarder sur",nj+nb-1,"soit ",(nb-1)/(nj+(nb-1))*100,"%")
