#!/usr/bin/python
#-*-coding:utf-8-*-

import random
from Bio.Seq import Seq
from statistics import *
import numpy as np
import sys
from Bio import SeqIO
import os
from Bio.SeqRecord import SeqRecord

def verifie(x):
	verif=open(x,'r')
	ligne=verif.readlines()
	verif.close()
	if len(ligne)==1:
		return ligne


for i in sys.argv:
	if i =="-out":
		out=sys.argv[sys.argv.index(i)+1]
	if i =="-data":
		data=sys.argv[sys.argv.index(i)+1]
	if i =="-draft":
		draft=sys.argv[sys.argv.index(i)+1]

tempo_out=out+"tempo/"

for rec in SeqIO.parse(data, "fasta"):
	if rec.name =="5":
		SeqIO.write(rec,tempo_out+"d5.fasta","fasta")
	if rec.name =="3":
		SeqIO.write(rec,tempo_out+"d3.fasta","fasta")


taille_data=len(SeqIO.read(tempo_out+"d5.fasta", "fasta"))
if taille_data <100 :
	type_blast="-task blastn-short"
else:
	type_blast=""
#-evalue 0.001 -qcov_hsp_perc 90 -perc_identity 90

os.system("""makeblastdb -in """+draft+""" -out """+tempo_out+"""draft -dbtype 'nucl'""")
#os.system("""blastn -query """+tempo_out+"""d5.fasta -db """+tempo_out+"""draft -out """+tempo_out+"""data.csv """+type_blast+""" -outfmt "10 sseqid sstrand" -evalue 0.00001""")
os.system("""blastn -query """+tempo_out+"""d5.fasta -db """+tempo_out+"""draft -out """+out+"""sens_5.html """+type_blast+""" -html -evalue 0.00001""")
# sens_5=verifie(tempo_out+"data.csv")


os.system("""makeblastdb -in """+draft+""" -out """+tempo_out+"""draft -dbtype 'nucl'""")
#os.system("""blastn -query """+tempo_out+"""d3.fasta -db """+tempo_out+"""draft -out """+tempo_out+"""data.csv """+type_blast+""" -outfmt "10 sseqid sstrand" -evalue 0.00001""")
os.system("""blastn -query """+tempo_out+"""d3.fasta -db """+tempo_out+"""draft -out """+out+"""sens_3.html """+type_blast+""" -html -evalue 0.00001""")
# sens_3=verifie(tempo_out+"data.csv")

# for rec in SeqIO.parse(data, "fasta"):
# 	if sens_3:
# 		if rec.name==sens_3[0]:
# 			if sens_3[1]=="plus"
# 				data_3=SeqRecord(Seq(rec.seq))
# 				data_3.id("3")
# 			if sens_3[1]=="minus"
# 				data_3=SeqRecord(Seq(rec.seq).reverse_complement())
# 				data_3.id("3")
