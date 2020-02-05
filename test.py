#!/usr/bin/python
#-*-coding:utf-8-*-
import os
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from statistics import median
import sys
import platform
import random


draft="/Users/Florent/pipeline/donnee_pour_git/result/assembly.fasta"
tempo_out="/Users/Florent/pipeline/donnee_pour_git/result/"
for rec in SeqIO.parse(draft, "fasta"):
	if rec.name == "5":
		SeqIO.write(rec,tempo_out+"1f.fasta","fasta")
		SeqIO.write(rec.reverse_complement(),tempo_out+"1fr.fasta","fasta")
	if rec.name =="3":
		SeqIO.write(rec,tempo_out+"2f.fasta","fasta")
		SeqIO.write(rec.reverse_complement(),tempo_out+"2fr.fasta","fasta")
