# LightMetaEbwt 

LightMetaEbwt is a new lightweight alignment-free and assembly-free framework for metagenomic classification that is combinatorial by nature and allows us to use little internal memory. Let *S* be a large collection of biological sequences comprising both reads and genomes. For simplicity, we denote by *R* the subset of reads (metagenomic sample) and by *G* the subset of genomes (reference database).

It takes in input:
- the extended Burrows–Wheeler transform (ebwt), or multi-string BWT, of collection *S*;
- the longest common prefix array (lcp) of collection *S*;
- the document array (da) of collection *S*.

By establishing a similarity degree between sequences that exploits the underlying properties of the eBWT, LightMetaEbwt performs a metagenomic classification by assigning any read in *R* to a unique genome in *G*. The method steps can be summarized in three steps, as follows: 

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

There are mainly two options to obtain the required data structures (ebwt, lcp, da) for *S*:
- one could build ebwt, lcp, and da starting from a fasta file containing both the reads in *R* and the genomes in *G*;
- one could build the data structures ebwt, lcp, and da separately for the set *G* and the set *R*, and then merge them to obtain ebwt, lcp, and da for the collection *S*.

The advantage of the latter choice lies in the fact that we can build the data structures of *G* only once if the set *G* of genomes is the same for each experiment.

To build ebwt, lcp, and da files from scratch from a single fasta file, one could use BCR [https://github.com/giovannarosone/BCR_LCP_GSA], or egsa [https://github.com/felipelouza/egsa] for instance. Note that egsa tool returns the three datastructures in a single file (fastaFile.K.gesa, with K being the number of sequences). The executable file EGSAtoBCR is to convert fastaFile.K.gesa into fileFasta.ebwt, fileFasta.lcp, and fileFasta.da -- use command EGSAtoBCR filefasta K.

To merge the data structures ebwt and lcp associated with *R* and *G*, one could use eGap [https://github.com/felipelouza/egap] and set the option -d to obtain the document array (da) of the merge. Note that in this case the output document array contains only 0s and 1s, and as to obtain the da of the entire collection *S* it is necessary to replace 0s with the values in da(*R*) and 1s with da(*G*).

On the other hand, exploiting the mathematical properties of the permutation associated with the
ebwt and lcp, one could use BCR [https://github.com/giovannarosone/BCR_LCP_GSA] incrementally in order to update the data structures obtained for *G* and obtain the data structures for the entire collection *S* (without constructing the eBWT and lcp for *R* from scratch) .

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

In Datasets, we provide some examples of simulated metagenomic samples. (See for details Datasets/Experiments_links.txt).

The dataset setB2 is a sample of 20,249,373 paired end short reads (100 bps) stored in setB2_1.fasta and setB2_2.fasta.

As preprocessing, we construct the datastructures ebwt, lcp, and da for the set *G* -- Refs.fasta is the fasta file of reference genomes (number of genomes: 930). 
Then, we merge the datastructures (ebwt, lcp, da) associated with Refs.fasta to those associated with the sets of reads (setB2_1.fasta and setB2_2.fasta, and setB2_1_RC.fasta and setB2_2_RC.fasta, which contain the reverse complement of setB2_1.fasta and setB2_2.fasta) as to obtain the datastructures for the four collections: 
setB2_1+Refs.fasta, setB2_1_RC+Refs.fasta, setB2_2+Refs.fasta, setB2_2_RC+Refs.fasta.

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
## References

[1] V. Guerrini and G. Rosins. Lightweight Metagenomic Classification via eBWT. Alcob 2019. LNCS, vol 11488, pp 112-124.

## Thanks

<small> Supported by the project Italian MIUR-SIR [CMACBioSeq][240fb5f5] ("_Combinatorial methods for analysis and compression of biological sequences_") grant n.~RBSI146R5L. P.I. Giovanna Rosone</small>

[240fb5f5]: http://pages.di.unipi.it/rosone/CMACBioSeq.html
