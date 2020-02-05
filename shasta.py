#!/usr/bin/python
#-*-coding:utf-8-*-
import os
from Bio.Seq import Seq
from Bio import SeqIO
from statistics import median
import sys
import platform

# reads="/Users/Florent/pipeline/reads/"

for i in sys.argv:
	if i =="-out":
		out=sys.argv[sys.argv.index(i)+1]
	if i =="-reads":
		reads=sys.argv[sys.argv.index(i)+1]
	if i =="-amorce":
		amorce=sys.argv[sys.argv.index(i)+1]
	if i =="-fichier":
		fichier=sys.argv[sys.argv.index(i)+1]
	if i =="-len":
		len2=sys.argv[sys.argv.index(i)+1]


pipenv=os.getcwd()
tempo=reads+"tempo/"
tempo_out=out+"tempo/"
tailleread=int(len2)

system=platform.system()
if system=="Linux":
	shasta_os="shasta-Linux-0.1.0"
if system=="Darwin":
	shasta_os="shasta-macOS-0.1.0"

premier='python3 fastaQtofasta.py -reads '+reads+" -out "+out+" -len "+len2
os.system(premier)
os.system("""makeblastdb -in """+tempo_out+"""multifasta.fasta -out """+tempo_out+"""target -dbtype 'nucl'""")
os.system("""blastn -query """+amorce+""" -db """+tempo_out+"""target -out """+tempo_out+"""read.fasta -task blastn-short -outfmt "10 sstrand sseqid" -evalue 0.001 -max_target_seqs 999999999""")
os.system("python3 multifastaPosttarget.py"+" -out "+out)
os.system("python3 selection_Q.py"+" -out "+out)
#os.system("python3 data.py")
trim=9
taille_contig=str()
ref_ok=False
cible_shasta=tempo_out+"read_taille.fasta"
un=tempo_out+"1.fasta"
deux=tempo_out+"2.fasta"



#while ref_ok is False :
while trim <19 and ref_ok is False:
	try:
		longueur_ref=sum([len(x.seq) for x in SeqIO.parse(tempo_out+"1f.fasta", "fasta")]+[len(x.seq) for x in SeqIO.parse(tempo_out+"2f.fasta", "fasta")])
	except:
		longueur_ref=sum([len(x.seq) for x in SeqIO.parse(un, "fasta")]+[len(x.seq) for x in SeqIO.parse(deux, "fasta")])
	# ref=open("ref.csv","a")
	# ref.write(str(longueur_ref)+"\n")
	# ref.close()
	shasta=pipenv+"/"+shasta_os+" --input "+cible_shasta+" --Reads.minReadLength "+str(tailleread)+" --Align.maxTrim "+str(trim)+" --output "+tempo_out+"Shastarun"
	os.system(shasta)
	mv="mv "+tempo_out+"Shastarun/Assembly.fasta "+tempo_out
	os.system(mv)
	rm="rm"+" -r "+tempo_out+"Shastarun"
	os.system(rm)
	assembly=tempo_out+"Assembly.fasta"
	nbcontig=len([x for x in SeqIO.parse(assembly, "fasta")])
	longueur_A=sum([len(x.seq) for x in SeqIO.parse(assembly, "fasta")])
	if nbcontig==2 and longueur_A>longueur_ref:
		num_contig=1
		for rec in SeqIO.parse(assembly, "fasta"):
			nom=tempo_out+str(num_contig)+'.fasta'
			SeqIO.write(rec,nom,"fasta")
			num_contig=num_contig+1
		os.system("""makeblastdb -in """+un+""" -out """+tempo_out+"""target -dbtype 'nucl'""")
		os.system("""blastn -query """+deux+""" -db """+tempo_out+"""target -out """+tempo_out+"""contig.fasta  -outfmt "10 sstrand sseqid"  -evalue 0.01""")
		tailletest=open(tempo_out+"contig.fasta","r")
		tailletest2=tailletest.readlines()
		tailletest.close()
		trim=trim+1
		if len(tailletest2) ==0:
			ref_ok=True
			longueur_ref=longueur_A
			cp1="cp "+un+" "+tempo_out+"1f.fasta"
			os.system(cp1)
			cp2="cp "+deux+" "+tempo_out+"2f.fasta"
			os.system(cp2)
			fichier_a=SeqIO.read(un, "fasta")
			taille_contig=str(len(fichier_a.seq))+","
			fichier_a=SeqIO.read(deux, "fasta")
			taille_contig=taille_contig+str(len(fichier_a.seq))
			robustesse=open(tempo_out+"robustesse.csv","a")
			ligne=str(ref_ok)+","+fichier+","+str(trim)+","+str(longueur_ref)+","+taille_contig+"\n"
			robustesse.write(ligne)
			robustesse.close()
	else :
		trim = trim +1
	#if trim ==20 and cible_shasta=="read_taille.fasta":
		#trim=9
	#	cible_shasta="readsurcible.fasta"
	#	trim=9
	#if trim ==20 and cible_shasta=="readsurcible.fasta":
	#	ref_ok=True
	#	ref=False



#ligne=str(rep)+","+str(ref)+","+str(ref_ok)+","+taille_contig+"\n"

# robustesse=open(tempo_out+"robustesse.csv","a")
# ligne=str(ref_ok)+","+fichier+","+str(longueur_ref)+","+taille_contig+"\n"
# robustesse.write(ligne)
# robustesse.close()
