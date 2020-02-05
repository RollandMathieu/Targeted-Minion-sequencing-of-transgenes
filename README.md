# Targeted-Minion-sequencing-of-transgenes
This pipeline allows de novo assembly of sequences flanking an amplified sequence according to the technique described in Boutigny et al. "Targeted Minion sequencing of transgenes".

It works only with MacOS/Linux.

### Installing
Installation needs 

conda https://docs.conda.io/en/latest/miniconda.html. Download python version 3.7 

blast ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.8.1/ for Mac and Linux sudo apt-get install ncbi-blast+
```
Conda create -n medaka -c bioconda -c conda-forge medaka 
pip install biopython 
git clone https://github.com/RollandMathieu/Targeted-Minion-sequencing-of-transgenes
cd Targeted-Minion-sequencing-of-transgenes 
```
Version linux : chmod ugo+x shasta-Linux-0.1.0

Version Mac : chmod ugo+x shasta-macOS-0.1.0

### Running
```
python3 shasta_run.py -len minimum_size_of_the_sequence -reads folder_containing_the_fastq_files -out result_folder -amorce file_containing_the_primers
```
A shasta.fasta file containing the 2 contigs is created. User has to determine which contig is at the 3’ and which contig is at the 5’ end of the primers. To do so, user has 3 possibilities : 
* blast, 
* perform a sanger sequencing on the post PCR product with one of the 2 primer,
* use the script blast.py
```
python3 blast.py -out result_folder -data file_containing_the_maximum_information (exemple sequence of the promotor) -draft file_shasta.py_generated_with_shasta_run.py
```
See the file data.fasta available as an exemple

The following part allows to complete the contig and to realize the polishing with medaka 
```
python3 pipeline.py -len 3000 -reads folder_containing your_files_fastq -out result_folder –amorce file_containing_your_primers -draft file_shasta.py_generated_with_shasta_run.py
```
In the file_shasta.py_generated_with_shasta_run.py, the orientation of the two contigs has to be indacted as a name >5 and >3 (see example shasta_orientated.fasta)

The draft assembly can be found in the folder out under the name draft.fasta, reads for the next step are in the same folder under the name read.fasta

To run the pollishing
```
conda activate medaka
medaka_consensus -i path_to_the_reads -d path_to_the_draft -o result_folder
```
The file consensus.fasta corresponds to the final result the expected quality is around Q30

### Authors
**Florent Fioriti** - *Initial work*
