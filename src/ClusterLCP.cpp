#include "Tools.h"

/*
ClusterLCP DETECTS CLUSTERS OF eBWT BY USING LCP AND DA VALUES

As preprocessing, it need to have the two data structures fileFasta.lcp and fileFasta.da computed.

Input: fileFasta, total number of reads, total number of genomes, and minimum value for LCP (alpha)
		
Output: a file with extension .clrs containing pairs ElementCluster (pStart,len) for each alpha-cluster detected

*/

dataTypeNChar clusterLCP(FILE* InFileLCP, FILE* InFileDA, FILE* fdCluster, dataTypelenSeq minimumLCP, dataTypeNSeq nRead, dataTypeNChar &maxLen)
{
	dataTypeNChar numClust=0;
	
    //Output file contains a collection of pairs (pStart, len) each one corresponding to a alpha-cluster
    ElementCluster cluster;
	cluster.pStart=0;
	cluster.len=0;
    
    //Read LCP and DA files
    dataTypeNChar numcharLCP;
    dataTypelenSeq* bufferLCP = new dataTypelenSeq[BUFFERLCPSIZE];
    
	dataTypeNChar numcharDA;
    dataTypeNSeq *bufferEle = new dataTypeNSeq[BUFFERLCPSIZE+1];
	
	fseek(InFileLCP, 0, SEEK_SET);
	fseek(InFileDA, 0, SEEK_SET);
    
    numcharLCP = fread(bufferLCP,sizeof(dataTypelenSeq),BUFFERLCPSIZE,InFileLCP);
    
	bufferEle[0]=0;
    numcharDA = fread(bufferEle+1,sizeof(dataTypeNSeq),BUFFERLCPSIZE,InFileDA);
    
    assert(numcharLCP==numcharDA);
	
    dataTypeNChar index=0; //position in the list of sorted suffixes
    
	bool init=false; //init is true if a cluster is open
	bool nR=false;	//nR is true when at least a read is in
	bool nG=false;	//nG is true when at least a genome is in
    
	//Start Clustering
	while ( (numcharLCP >0 ) && ( numcharDA >0 )  ) {
        
        for(dataTypeNChar indexbuffer=0; indexbuffer<numcharLCP; indexbuffer++)
        {
			if(bufferLCP[indexbuffer]>=minimumLCP) //Start or remain in a cluster
            {
                if (not init) //Start a new cluster
                {
					init=true;
                    cluster.pStart = index-1;
					if (bufferEle[indexbuffer]<nRead)
						nR=true;
					else
						nG=true;
                }
				
				if (bufferEle[indexbuffer+1]<nRead)
					nR=true;
				else
					nG=true;
			}
            else    //End a cluster
			{
				if (init && nR && nG)
				{
                    cluster.len = index - cluster.pStart;
                    if (maxLen<cluster.len)
                        maxLen=cluster.len;

					fwrite(&cluster,sizeof(ElementCluster),1,fdCluster);
                    numClust++;
				}
                
				init=false;
				nR=false;
				nG=false;
				
				cluster.pStart=0;
				cluster.len=0;
			}
			index++;
        }  //end-for
		
        //Read the LCP and DA files
        bufferEle[0]=bufferEle[numcharDA];
        numcharLCP = fread(bufferLCP,sizeof(dataTypelenSeq),BUFFERLCPSIZE,InFileLCP);
        numcharDA = fread(bufferEle+1,sizeof(dataTypeNSeq),BUFFERLCPSIZE,InFileDA);
        assert(numcharLCP==numcharDA);
        
    }//end-while
    
    
    //We need to close a possible final cluster remained open iff it contains a read and a genome
    if(init && nR && nG)
    {
		cluster.len = index - cluster.pStart;
        if (maxLen<cluster.len)
            maxLen=cluster.len;

        fwrite(&cluster,sizeof(ElementCluster),1,fdCluster);
        numClust++;
    
    }
	delete[] bufferLCP;
	delete[] bufferEle;
    
    return numClust;
}



int main(int argc, char **argv) {

    time_t t_total=0;
    clock_t c_total=0;

	if( argc != 5 )
    {
      std::cerr << "Error usage: " << argv[0] << " fileFasta numReads numGenomes alpha" << std::endl;
       exit(1);
    }
    
    string fileFasta=argv[1]; 
	dataTypelenSeq alpha;
    dataTypeNSeq numReads, numGenomes;
	
    sscanf(argv[2], "%u", &numReads);
	sscanf(argv[3], "%u", &numGenomes);
    sscanf(argv[4], "%u", &alpha);
    
    string fnLCP, fnDA, fileOutput;
	fnLCP = fileFasta+".lcp\0";
    fnDA =fileFasta+".da\0";
    std::stringstream ssout;
    ssout << fileFasta << "." << alpha << ".clrs\0";
	fileOutput=ssout.str();
    
    FILE *OutCluster=fopen(fileOutput.c_str(), "w");
    if(OutCluster==NULL) {
        cerr << "Error opening " << fileOutput << ".";
        exit(1);
    }
	
	//Open LCP and DA files
	FILE * InLCP, *InDA;
    
    InLCP = fopen(fnLCP.c_str(), "rb");
    if (InLCP==NULL) {
        std::cerr << "Error opening " << fnLCP << "." << std::endl;
        printf("fopen failed, errno = %d\n", errno);
        exit (EXIT_FAILURE);
    }
    
    
    InDA = fopen(fnDA.c_str(), "rb");
    if (InDA==NULL) {
        std::cerr << "Error opening " << fnDA << "." << std::endl;
        printf("fopen failed, errno = %d\n", errno);
        exit (EXIT_FAILURE);
    }
    
	//Clustering
    time_start(&t_total, &c_total);
    dataTypeNChar maxLen=0;
    
    dataTypeNChar nClusters=clusterLCP(InLCP,InDA,OutCluster, alpha, numReads, maxLen);

    fclose(InLCP);
    fclose(InDA);
    fclose(OutCluster);
    
    string fileaux=fileFasta.substr(0,fileFasta.find(".fasta"))+".out";

	//Write auxiliary file
    FILE * outAux = fopen(fileaux.c_str(), "w");
    if (outAux==NULL) {
        std::cerr << "Error opening " << fileaux << "." << std::endl;
        printf("fopen failed, errno = %d\n", errno);
        exit (EXIT_FAILURE);
    }
    
	fwrite(&numReads,sizeof(dataTypeNSeq),1,outAux);
    fwrite(&numGenomes,sizeof(dataTypeNSeq),1,outAux);
    fwrite(&alpha,sizeof(dataTypelenSeq),1,outAux);
    fwrite(&maxLen,sizeof(dataTypeNChar),1,outAux);
                                    
    fclose(outAux);
                                    
    cout << "Clustering process with alpha=" << alpha << " completed.\nTotal number of clusters: " << nClusters << ".\nMaximum cluster size: " << maxLen << "." << endl;
	fprintf(stdout,"Time: %.6lf\n", time_stop(t_total, c_total));
	return 0;
}
