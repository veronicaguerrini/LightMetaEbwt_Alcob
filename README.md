# LightMetaEbwt 

LightMetaEbwt is a new lightweight alignment-free and assembly-free framework for metagenomic classification that is combinatorial by nature and allows us to use little internal memory. Let *S* be a large collection of biological sequences comprising both reads and genomes. For simplicity, we denote by *R* the subset of reads (metagenomic sample) and by *G* the subset of genomes (reference database).

It takes in input:
- the extended Burrows–Wheeler transform (ebwt), or multi-string BWT, of collection *S*;
- the longest common prefix array (lcp) of collection *S*;
- the document array (da) of collection *S*.

By establishing a similarity degree between sequences that exploits the underlying properties of the eBWT, LightMetaEbwt performs a metagenomic classification by assigning any read in *R* to a unique genome in *G*. The method steps can be summarised in three steps, as follows: 

(1) we detect and keep some blocks of ebwt(*S*) in which the associated suffixes share a common context of a minimum length *alpha*, and to which both reads and genomes belong; 

(2) we analyze these interesting blocks in order to evaluate a degree of similarity between any read and any genome in *S*, and we discard similarity values that do not exceed a threshold value *beta*; 

(3) we perform the read assignment: for every read, either we retrieve the unique genome of belonging or we report that it is not possible to identify it.


### Install

```sh
git clone https://github.com/veronicaguerrini/LightMetaEbwt
cd LightMetaEbwt
install.sh
```
### Preprocessing step

The required data structures can be obtained by running BCR [https://github.com/giovannarosone/BCR_LCP_GSA] on input the fasta file of the entire collection *S*.

### Run

The three steps are accomplished by running:

(1) ClusterLCP with input parameters name of the fasta file, total number of reads in *S*, total number of genomes in *S* and *alpha*;

(2) ClusterBWT with input parameters name of the fasta file, length of reads, and *beta*;

(3) Classify providing in input txt files obtained by running ClusterBWT, the total number of genomes in *S*, and name of the output file.

```sh
 ClusterLCP fileFasta numReads numGenomes alpha
 ClusterBWT fileFasta readLength beta
 Classify N fileInput1 ... fileInputN numGenomes fileOutput
```
Recall that in order to run ClusterLCP we need to have fileFasta.lcp and fileFasta.da computed, while to run ClusterBWT we need fileFasta.da and fileFasta.ebwt.

### Example

SetB2 is a set of 20,249,373 paired end short reads (100 bps).
We construct datastructures (ebwt, lcp, da) for setB2_1.fasta and setB2_2.fasta, and for setB2_1_RC.fasta and setB2_2_RC.fasta, which contain their reverse complement.
We construct datastructures (ebwt, lcp, da) for the set *G* -- Refs.fasta is the fasta file of reference genomes (number of genomes: 930). (See for details Datasets/Experiments_links.txt). 
We merge the datastructures (ebwt, lcp, da) associated with Refs.fasta to those associated with the sets of reads, as to obtain the datastructures for the four collections: setB2_1+Refs.fasta, setB2_1_RC+Refs.fasta, setB2_2+Refs.fasta, setB2_2_RC+Refs.fasta.

We assign any read (or its reverse complement) in setB2 to a reference genome in *G*.
For our analysis, we set the minimum common context alpha=16 and the minimum similarity score beta=0.25.

```sh
for X in 1 1_RC 2 2_RC
do
 ClusterLCP setB2_$X+Refs.fasta 20249373 930 16
 ClusterBWT setB2_$X+Refs.fasta 100 0.25 
done
```
The output of steps (1) and (2) are then given in input to step (3).

```sh
Classify 4 Clustering_results_B0.25_setB2_1+Refs.fasta.txt Clustering_results_B0.25_setB2_1_RC+Refs.fasta.txt Clustering_results_B0.25_setB2_2+Refs.fasta.txt Clustering_results_B0.25_setB2_2_RC+Refs.fasta.txt 930 results_setB2.txt
```
---
<small> Supported by the project Italian MIUR-SIR [CMACBioSeq][240fb5f5] ("_Combinatorial methods for analysis and compression of biological sequences_") grant n.~RBSI146R5L. P.I. Giovanna Rosone</small>

[240fb5f5]: http://pages.di.unipi.it/rosone/CMACBioSeq.html
