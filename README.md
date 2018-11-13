# LightMetaEbwt 

LightMetaEbwt is a new lightweight alignment-free and assembly-free framework for metagenomic classification that compares each unknown sequence in the sample to a collection of known genomes.

It takes in input:
- the extended Burrows–Wheeler transform (ebwt), or multi-string BWT,
- the longest common prefix array (lcp),
- the document array (da),

of a very large collection *S* of strings including both reads and genomes. The required data structures can be obtained by running BCR [https://github.com/giovannarosone/BCR_LCP_GSA] on input the fasta file of the entire collection *S*.

### Install

```sh
git clone https://github.com/veronicaguerrini/LightMetaEbwt
cd LightMetaEbwt
install.sh
```

### Run
Our method works in three steps: 

(1) we detect and keep some blocks of ebwt(*S*) in which the associated suffixes share a common context of a minimum length *alpha*, and to which both reads and genomes belong; 

(2) we analyze these interesting blocks in order to evaluate a degree of similarity between any read and any genome in *S*, and we discard similarity values that do not exceed a threshold value *beta*; 

(3) we perform the read assignment: for every read, either we retrieve the unique genome of belonging or we report that it is not possible to identify it.

The three steps are accomplished by running:

(1) ClusterLCP with input parameters name of the fasta file, total number of reads in *S*, total number of genomes in *S* and *alpha*;

(2) ClusterBWT with input parameters name of the fasta file, length of reads, and *beta*;

(3) Classify providing in input txt files obtained by running ClusterBWT, and the total number of genomes in *S*.

```sh
 ClusterLCP fileFasta numReads numGenomes alpha
 ClusterBWT fileFasta readLength beta
 Classify N fileInput1 ... fileInputN numGenomes
```
Note that in order to run ClusterLCP we need to have fileFasta.lcp and fileFasta.da computed, while to run ClusterBWT we need fileFasta.da and fileFasta.ebwt.

### Example
SetB2_1+Refs.fasta and setB2_2+Refs.fasta are two sets of sequences containing paired end reads (number of reads in each set: 20,249,373) and reference genomes (number of genomes: 930). (See for details Datasets/Experiments_links.txt). 

We set alpha=16 and beta=0.36, and then we classify each read comparing similarity values reported by both ends of the same read.

```sh
./ClusterLCP setB2_1.fasta 20249373 930 16
./ClusterBWT setB2_1.fasta 100 0.36
./ClusterLCP setB2_2.fasta 20249373 930 16
./ClusterBWT setB2_2.fasta 100 0.36
./Classify 2 Clustering_results_B0.36_setB2_1.fasta.txt Clustering_results_B0.36_setB2_2.fasta.txt 930
```
---
<small> Supported by the project Italian MIUR-SIR [CMACBioSeq][240fb5f5] ("_Combinatorial methods for analysis and compression of biological sequences_") grant n.~RBSI146R5L. P.I. Giovanna Rosone</small>

[240fb5f5]: http://pages.di.unipi.it/rosone/CMACBioSeq.html
