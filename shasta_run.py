#!/usr/bin/python
#-*-coding:utf-8-*-
import os
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from statistics import median
import sys
import random

reads="/Users/Florent/pipeline/reads/"
tempo="/Users/Florent/pipeline/reads/tempo/"

for i in sys.argv:
	if i =="-out":
		out=sys.argv[sys.argv.index(i)+1]
	if i =="-reads":
		reads=sys.argv[sys.argv.index(i)+1]
	if i =="-amorce":
		amorce=sys.argv[sys.argv.index(i)+1]
	if i =="-len":
		len2=sys.argv[sys.argv.index(i)+1]

os.system("mkdir "+out)
tempo=reads+"tempo/"
tempo_out=out+"tempo/"
os.system("mkdir "+tempo)
temporaire="mkdir "+tempo_out
os.system(temporaire)
robu=out+"tempo/robustesse.csv"
robustesse=open(robu,"w")
robustesse.close()
### obtenir cette liste automatiquement
# liste={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}
reset="mv "+reads+"*.fastq "+tempo
os.system(reset)
merge=os.listdir(tempo)


SeqIO.write(SeqRecord(Seq("A")), tempo_out+"1.fasta", "fasta")
SeqIO.write(SeqRecord(Seq("A")), tempo_out+"2.fasta", "fasta")
for nb in range(1,2):
	fastq=[x for x in merge if x[-5:] =="fastq"]
	for l in range(1,int(len(fastq)+1)):
		reset="mv "+reads+"*.fastq "+tempo
		os.system(reset)
		run=random.sample(fastq,nb)
		[fastq.remove(x) for x in run]
		for i in run:
			move="mv "+tempo+i+" "+reads
			os.system(move)
		shasta="python3 shasta.py -fichier "+"/".join(run)+" -out "+out+" -reads "+reads+" -amorce "+amorce+" -len "+len2
		os.system(shasta)

nb=len(fastq)
# reset="mv /Users/Florent/pipeline/reads/*.fastq /Users/Florent/pipeline/reads/tempo/"
os.system(reset)
run=random.sample(fastq,nb)
for i in run:
	move="mv "+tempo+i+" "+reads
	os.system(move)
shasta="python3 shasta.py -fichier "+"/".join(run)+" -out "+out+" -reads "+reads+" -amorce "+amorce+" -len "+len2
os.system(shasta)
reset="mv "+reads+"tempo/*.fastq "+reads
os.system(reset)

supp_tempo="rm -rf "+tempo[:-1]
os.system(supp_tempo)


os.system("cat "+tempo_out+"1f.fasta "+tempo_out+"2f.fasta >>"+out+"shasta.fasta")

# R1=SeqIO.read(tempo_out+"1f.fasta","fasta")
# R1r=R1.reverse_complement()
# SeqIO.write(R1r,tempo_out+"1fr.fasta","fasta")
#
# R2=SeqIO.read(tempo_out+"2f.fasta","fasta")
# R2r=R2.reverse_complement()
# SeqIO.write(R2r,tempo_out+"2fr.fasta","fasta")
