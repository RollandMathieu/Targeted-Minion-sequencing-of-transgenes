#!/usr/bin/python
#-*-coding:utf-8-*-

import os
import sys
from Bio.Seq import Seq
from Bio import SeqIO
import platform

def verifie(x):
	verif=open(x,'r')
	ligne=verif.readlines()
	verif.close()
	for i in ligne:
		ex=i.split(",")
		if ex[0]=="plus" and int(ex[3])==1 and int(ex[4]) < 1000 and int(ex[5]) > 1000 :
			return i

def verifieR(x):
	verif=open(x,'r')
	ligne=verif.readlines()
	verif.close()
	for i in ligne:
		ex=i.split(",")
		if ex[0]=="plus" and int(ex[3])==int(ex[4]) and int(ex[7])-int(ex[5]) < 1000 and int(ex[6]) >1000:
			return i

pipenv=os.getcwd()

system=platform.system()
if system=="Linux":
	chemin_spades=pipenv+"/spades_linux/bin/spades.py"
if system=="Darwin":
	chemin_spades=pipenv+"/spades_macos/bin/spades.py"


for i in sys.argv:
	if i =="-reads":
		reads=sys.argv[sys.argv.index(i)+1]
	if i =="-out":
		out=sys.argv[sys.argv.index(i)+1]
	if i =="-amorce":
		amorce=sys.argv[sys.argv.index(i)+1]
	if i =="-len":
		len2=sys.argv[sys.argv.index(i)+1]
	if i =="-draft":
		draft=sys.argv[sys.argv.index(i)+1]

tempo_out=out+"tempo/"
run="mkdir "+out
os.system(run)
temporaire="mkdir "+tempo_out
os.system(temporaire)

for rec in SeqIO.parse(amorce, "fasta"):
	if rec.name =="5":
		SeqIO.write(rec,tempo_out+"5.fasta","fasta")
	if rec.name =="3":
		SeqIO.write(rec,tempo_out+"3.fasta","fasta")

for rec in SeqIO.parse(draft, "fasta"):
	if rec.name == "3":
		SeqIO.write(rec,tempo_out+"1f.fasta","fasta")
		SeqIO.write(rec.reverse_complement(),tempo_out+"1fr.fasta","fasta")
	if rec.name =="5":
		SeqIO.write(rec,tempo_out+"2f.fasta","fasta")
		SeqIO.write(rec.reverse_complement(),tempo_out+"2fr.fasta","fasta")
###pour remplir le coter 3'

premier='python3 fastaQtofasta.py -reads '+reads+" -out "+out+" -len "+len2
os.system(premier)
os.system("""makeblastdb -in """+tempo_out+"""multifasta.fasta -out """+tempo_out+"""target -dbtype 'nucl'""")
os.system("""blastn -query """+tempo_out+"""3.fasta -db """+tempo_out+"""target -out """+tempo_out+"""read.fasta -task blastn-short -outfmt "10 sstrand sseqid" -evalue 0.001 -max_target_seqs 999999999""")
os.system("python3 multifastaPosttarget.py"+" -out "+out)
os.system("python3 selection_Q.py"+" -out "+out)

fait=int()
a=str()
#début de la boucle ici
while not a:
	os.system("python3 data.py"+" -out "+out)
	os.system("python3 splitfasta.py"+" -out "+out)

	chemin_spades_source=tempo_out+"/readsplit.fasta"

	spades="python3 "+chemin_spades+" -s "+chemin_spades_source+" -o "+tempo_out+"assembly --only-assembler -m 4 -t 2"
	os.system(spades)

	scafold="mv "+tempo_out+"/assembly/scaffolds.fasta "+tempo_out
	os.system(scafold)

	os.system("""makeblastdb -in """+tempo_out+"""scaffolds.fasta -out """+tempo_out+"""target -dbtype 'nucl'""")
	os.system("""blastn -query """+tempo_out+"""3.fasta -db """+tempo_out+"""target -out """+tempo_out+"""read.fasta -task blastn-short -outfmt "10 sstrand sseqid" -evalue 0.001 -max_target_seqs 999999999 -qcov_hsp_perc 90 -perc_identity 90""")
	os.system("python3 scaffoldstofasta.py"+" -out "+out)
	os.system("python3 contig.py"+" -out "+out)
	os.system("makeblastdb -in """+tempo_out+"""contigsurcible.fasta -out """+tempo_out+"""target -dbtype 'nucl'""")
	os.system("""blastn -query """+tempo_out+"""1f.fasta -db """+tempo_out+"""target -out """+tempo_out+"""assembly.csv  -outfmt "10 sstrand sseqid qstart sstart length" -evalue 0.0000001 -max_target_seqs 999999999""")
	a=verifie(tempo_out+"assembly.csv")
	os.system("""blastn -query """+tempo_out+"""1f.fasta -db """+tempo_out+"""target -out """+tempo_out+"""assembly.html  -html -evalue 0.0001 -max_target_seqs 999999999""")

	if a:
		fait=tempo_out+"1f.fasta"
	if not a :
		os.system("""blastn -query """+tempo_out+"""1fr.fasta -db """+tempo_out+"""target -out """+tempo_out+"""assembly.csv  -outfmt "10 sstrand sseqid qstart sstart length" -evalue 0.0000001 -max_target_seqs 999999999""")
		a=verifie(tempo_out+"assembly.csv")
		os.system("""blastn -query """+tempo_out+"""1fr.fasta -db """+tempo_out+"""target -out """+tempo_out+"""assembly.html  -html -evalue 0.0001 -max_target_seqs 999999999""")

		if a :
			fait=tempo_out+"1fr.fasta"

#fin de la boucle ici

#Partie qui récupère la partie manquante
debutp=open(tempo_out+"contigsurcible.fasta",'r')
ex=a.split(",")
contig=debutp.readlines()
compteur=0
while compteur < len(contig):
	if contig[compteur].startswith(">"):
		ligne2=contig[compteur].split(",")
		ligne=ligne2[1]
		sens=ligne2[0]
		ligne3=ligne2[0]
		if ligne[:-1]==ex[2] and ligne3[1:]==ex[1]:
			seq=contig[compteur+1]
			start=seq[:int(ex[4])-3]
	compteur=compteur+1

debutp.close()
resultat=open(tempo_out+"sequence.fasta","w")
resultat.write(">3'\n")
source=fait

for rec in SeqIO.parse(source, "fasta"):
	final1=start+str(rec.seq)
	resultat.write(final1)
resultat.write("\n")
resultat.close()
##################partie en 5'
os.system("""makeblastdb -in """+tempo_out+"""multifasta.fasta -out """+tempo_out+"""target -dbtype 'nucl'""")
os.system("""blastn -query """+tempo_out+"""5.fasta -db """+tempo_out+"""target -out """+tempo_out+"""read.fasta -task blastn-short -outfmt "10 sstrand sseqid" -evalue 0.001 -max_target_seqs 999999999""")
os.system("python3 multifastaPosttarget.py"+" -out "+out)
os.system("python3 selection_Q.py"+" -out "+out)


r5a=tempo_out+"2f.fasta"
r5b=tempo_out+"2fr.fasta"


fait=int()
a=str()
#début de la boucle ici
while not a:
	os.system("python3 data.py"+" -out "+out)
	os.system("python3 splitfasta.py"+" -out "+out)


	chemin_spades_source=tempo_out+"/readsplit.fasta"

	spades="python3 "+chemin_spades+" -s "+chemin_spades_source+" -o "+tempo_out+"assembly --only-assembler -m 4 -t 2"
	os.system(spades)

	scafold="mv "+tempo_out+"/assembly/scaffolds.fasta "+tempo_out
	os.system(scafold)

	os.system("makeblastdb -in "+tempo_out+"scaffolds.fasta -out "+tempo_out+"target -dbtype 'nucl'")
	os.system("""blastn -query """+tempo_out+"""5.fasta -db """+tempo_out+"""target -out """+tempo_out+"""read.fasta -task blastn-short -outfmt "10 sstrand sseqid" -evalue 0.001 -max_target_seqs 999999999 -qcov_hsp_perc 90 -perc_identity 90""")
	os.system("python3 scaffoldstofasta.py"+" -out "+out)
	os.system("python3 contig.py"+" -out "+out)
	os.system("makeblastdb -in "+tempo_out+"contigsurcible.fasta -out "+tempo_out+"target -dbtype 'nucl'")
	os.system("""blastn -query """+r5a+""" -db """+tempo_out+"""target -out """+tempo_out+"""assembly2.csv  -outfmt "10 sstrand sseqid qend qlen send length slen" -evalue 0.0000001 -max_target_seqs 999999999""")
	a=verifieR(tempo_out+"assembly2.csv")
	if a:
		fait=r5a
	if not a :
		os.system("""blastn -query """+r5b+""" -db """+tempo_out+"""target -out """+tempo_out+"""assembly2.csv  -outfmt "10 sstrand sseqid qend qlen send length slen" -evalue 0.0000001 -max_target_seqs 999999999""")## bon blast
		a=verifieR(tempo_out+"assembly2.csv")
		if a :
			fait=r5b

#fin de la boucle ici
debutp=open(tempo_out+"contigsurcible.fasta",'r')
ex=a.split(",")
contig=debutp.readlines()
compteur=0
while compteur < len(contig):
	if contig[compteur].startswith(">"):
		ligne2=contig[compteur].split(",")
		ligne=ligne2[1]
		ligne3=ligne2[0]
		sens=ligne2[0]
		if ligne[:-1]==ex[2] and ligne3[1:]==ex[1]:#and sens[1:]=="plus":
			seq=contig[compteur+1]
			fin=int(ex[7])-int(ex[5])
			start=seq[-fin:]
	compteur=compteur+1

debutp.close()
resultat=open(tempo_out+"sequence.fasta","a")
resultat.write(">5'\n")
source=fait

for rec in SeqIO.parse(source, "fasta"):
	final1=str(rec.seq)+start
	resultat.write(final1)
resultat.close()
premier='python3 fastaQtofasta.py -reads '+reads+" -out "+out+" -len "+len2
os.system(premier)

# os.system("conda activate medaka")
# medaka="medaka_consensus -i "+tempo_out+"multifasta.fasta -d "+tempo_out+"sequence.fasta -o "+out
# os.system(medaka)

move="mv "+tempo_out+"sequence.fasta "+out+"draft.fasta"
os.system(move)
move="mv "+tempo_out+"multifasta.fasta "+out+"reads.fasta"
os.system(move)
