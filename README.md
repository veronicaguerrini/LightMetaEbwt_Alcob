# LightMetaEbwt 

LightMetaEbwt is a new lightweight alignment-free and assembly-free framework for metagenomic classification that compares each unknown sequence in the sample to a collection of known genomes.

Takes in input:
- the (extended) Burrows–Wheeler transform (multi-string BWT)
- the longest common prefix array 
- the document array 

of a very large collection *S* of strings including both reads and genomes. The required data structures can be obtained by running BCR [https://github.com/giovannarosone/BCR_LCP_GSA] on input the fasta file of the entire collection *S*.

Soon online

### Install

```sh
git clone https://github.com/veronicaguerrini/LightMetaEbwt
cd LightMetaEbwt-master
```

```sh
./install.sh
```


### Run

```sh
 ./ClusterLCP fileFasta numReads numGenomes alpha
 ./ClusterBWT fileFasta readLength beta
 ./Classify numFile fileInput numGenomes
```

### Example
```sh
./ClusterLCP setB2_1.fasta 20249373 930 16
./ClusterLCP setB2_2.fasta 20249373 930 16
./ClusterBWT setB2_1.fasta 100 0.36
./ClusterBWT setB2_2.fasta 100 0.36
./Classify 2 Clustering_results_B0.36_setB2_1.fasta.txt Clustering_results_B0.36_setB2_2.fasta.txt 930
```


---
<small> Supported by the project Italian MIUR-SIR [CMACBioSeq][240fb5f5] ("_Combinatorial methods for analysis and compression of biological sequences_") grant n.~RBSI146R5L. P.I. Giovanna Rosone</small>

[240fb5f5]: http://pages.di.unipi.it/rosone/CMACBioSeq.html
