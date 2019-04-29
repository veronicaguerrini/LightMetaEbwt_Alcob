#include "Tools.h"

/*
ClusterBWT REFINES CLUSTERS OF eBWT BY USING BWT SYMBOLS

As preprocessing, it need to have the two data structures fileFasta.da and fileFasta.ebwt computed.

Input: fileFasta, read length, minimum similarity value beta
 
Output: A .txt file named Clustering_results_B(beta)_(fileFasta).txt having as many lines as the number of reads

Each line of the output file corresponds to a unique read R_j, and according to the minimum similarity value beta it can be either empty or storing the following information: 
maximum similarity value for R_j, the set S={(R_k, positive similarity between R_j and R_k), being R_k a reference genome}.
*/


dataTypeSim FindMinimumKeepUnmatched(dataTypeSim v1,dataTypeSim v2, dataTypeSim *v3, uint *diff_ref){
	dataTypeSim w=0;
	if (v1<=v2)
	{
		w=v1;
		*diff_ref+=v2-v1;
	}
	else //v1>v2 we use some placeholders
	{
		if (v1<=v2+*v3)
    {
			w=v1;
			*v3-=(v1-v2);
		}
		else
		{
			w=v2+*v3;   //v2+v3<v1<255
			*v3=0; //all placeholders used
		}
	}
	return w;
}

dataTypeSim FindMinimum(dataTypeSim v1,dataTypeSim v2, dataTypeSim *v3){
	//Case: #A(or #C,#G,#T) read v1 <= v2 #A(or #C,#G,#T) reference
	dataTypeSim w=0;
	if (v1<=v2)
		w=v1;
	else //v1>v2 we use some placeholders
	{
		if (v1<=v2+*v3)
		{
			w=v1;
			*v3-=(v1-v2);
		}
		else
		{
			w=v2+*v3;   //v2+v3<v1<255
			*v3=0; //all placeholders used
		}
	}
	return w;
}

dataTypeSim HowManyJolly(dataTypeSim v1, uint v2){
	if (v1<=v2)
		return v1;
	else
		return v2;
}



dataTypeNChar clusterRefine(FILE* InFileCluster, FILE* InFileDA, FILE* InFileBWT, dataTypeNSeq nRead, dataTypeNSeq tot, vector< vector<dataTypeSim> > &SimArray_)
{
    dataTypeNChar AnaCluster;
    
    //sigma is to store the BWT symbols
	dataTypedimAlpha *sigma;
  	sigma= new dataTypedimAlpha[USim_MAX];
  	dataTypedimAlpha index_letter_0=1;
  	dataTypedimAlpha index_letter_1=5;
  	for (uchar i=0; i<USim_MAX; i++)
 	{
		if( i==65 || i==67 || i==71 || i==84 )//A,C,G,T
		{
		sigma[i]=index_letter_0;
		index_letter_0++;
		}
		else if( i==66 || i==68 || i==72 || i==75 || i==77 || i==82 || i==83 || i==86 || i==87 || i==89 )//B,D,H,K,M,R,S,V,W,Y
		{
		sigma[i]=index_letter_1;
		index_letter_1++;
		}
		else //N,$
			sigma[i]=0;
	}
	
   dataTypeNSeq *mapID;
	mapID = new dataTypeNSeq[tot];
   dataTypeNSeq *mapIDinv;
	mapIDinv = new dataTypeNSeq[sizeMaxBuf];
    
	//To read InFileCluster
    dataTypeNChar numcharCluster;
	ElementCluster* clusterbuffer= new ElementCluster[BUFFERCLUSTER];
    dataTypeNSeq sizec, sizer, counter;
        
    //To read InFileDA
	dataTypeNSeq *elebuffer = new dataTypeNSeq[sizeMaxBuf];
    
	//To read InFileBWT
	dataTypedimAlpha *BWTbuffer= new dataTypedimAlpha[sizeMaxBuf];
	
    //cout << "Allocating memory..." << endl;
    
    dataTypeSim **CheckFreq; //Store the number of (A,C,G,T,placeholders) for each sequence
	CheckFreq=new dataTypeSim*[sizeMaxBuf];
    for (uint i = 0; i <sizeMaxBuf; ++i)
    	CheckFreq[i] = new dataTypeSim[index_letter_1];
    
    uint symbol;
    dataTypeNSeq entry;
    std::pair<std::set<dataTypeNSeq>::iterator,bool> ret;
    set<dataTypeNSeq> idInCluster;
    
    dataTypeSim t=0; //counter (t cannot be greater than 255, since each read has length 101)
    dataTypeSim placeholder=0; //to decrease placeholder symbols used
    uint diff_ref=0;//counter storing reference BWTletters that could match read placeholders
    
    //Read InFileCluster
	numcharCluster=fread(clusterbuffer,sizeof(ElementCluster),BUFFERCLUSTER,InFileCluster);
    
	AnaCluster=0;
    
    while (numcharCluster>0)
    {
        for (dataTypeNChar c=0;c<numcharCluster; c++)
        {
            
            //Read the BWT letters and DA file
            fseek(InFileBWT, (clusterbuffer[c].pStart)*sizeof(dataTypedimAlpha), SEEK_SET);
            fseek(InFileDA, (clusterbuffer[c].pStart)*sizeof(dataTypeNSeq), SEEK_SET);
            
            fread(BWTbuffer,sizeof(dataTypedimAlpha),clusterbuffer[c].len, InFileBWT);
            fread(elebuffer,sizeof(dataTypeNSeq),clusterbuffer[c].len, InFileDA);
            
            idInCluster.clear();
            sizer=0;
            sizec=0;
            
            for(dataTypeNChar i=0; i<clusterbuffer[c].len; i++)
            {
                ret=idInCluster.insert(elebuffer[i]);
				if ( (ret.second==true) && (elebuffer[i]<nRead) )
					sizer++;
            }
            sizec=idInCluster.size();
            
			assert(sizec>sizer);
            
			//cout << "We define maps..." << endl;
            
            counter=0;
            for (std::set<dataTypeNSeq>::iterator it=idInCluster.begin(); it!=idInCluster.end(); it++) //Define the map for sequences in cluster
            {
                mapID[*it]=counter;
                mapIDinv[counter]=*it;
               
                counter++;
            }
            assert(counter==sizec);
            idInCluster.clear();
            
            //Reset CheckFreq
            for(dataTypeNSeq i=0; i<sizec; i++)
            {
                for(uint k=0; k<index_letter_1; k++)
                    CheckFreq[i][k]=0;
            }
            
            
            for (dataTypeNChar i=0; i<clusterbuffer[c].len; i++)
            { 
                entry=mapID[elebuffer[i]];
                if (CheckFreq[entry][sigma[BWTbuffer[i]]]<USim_MAX)
              		CheckFreq[entry][sigma[BWTbuffer[i]]]++;
            }

            //cout << "Comparing..." << endl;
            
            for (dataTypeNSeq i=0; i<sizer;i++) //range over reads
            {
                for (dataTypeNSeq j=sizer; j<sizec;j++)//range over references
                {
			t=0;
                    diff_ref=0;
                    placeholder=CheckFreq[j][0];
			
                    for (dataTypedimAlpha k=1;k<index_letter_0; k++)
              		t+=FindMinimum(CheckFreq[i][k],CheckFreq[j][k], &placeholder); //Case A,C,G,T --> A,C,G,T,($,N)
		for (dataTypedimAlpha k=index_letter_0;k<index_letter_1; k++) //Case B,D,H,K,M,R,S,V,W,Y --> B,D,H,K,M,R,S,V,W,Y,($,N)
			t+=FindMinimumKeepUnmatched(CheckFreq[i][k],CheckFreq[j][k], &placeholder, &diff_ref);
			
           	 t+=HowManyJolly(CheckFreq[i][0],placeholder+diff_ref);//Case (N,$) --> iupac,$,N
					
                    uint temp = SimArray_[mapIDinv[i]][mapIDinv[j]-nRead]+t;
                    assert(temp <USim_MAX);
                    SimArray_[mapIDinv[i]][mapIDinv[j]-nRead]=temp;
                }
            }
            AnaCluster++;
            
        }//end-for
        
        numcharCluster=fread(clusterbuffer,sizeof(ElementCluster),BUFFERCLUSTER,InFileCluster);
    }//end-while
    
    delete[] clusterbuffer;
	delete[] elebuffer;
	delete[] BWTbuffer;
	
    delete[] mapID;
    delete[] mapIDinv;
	
	for(uint k= 0; k < sizeMaxBuf; ++k)
   	 	delete [] CheckFreq[k];
    
    return AnaCluster;
}

void clusterChoose(FILE* OutFileClusters, dataTypelenSeq &norm, float &beta, dataTypeNSeq nRead, dataTypeNSeq nRef, vector< vector<dataTypeSim> > &SimArray_)
{
    dataTypeSim maxSim=0; //maximum similarity value

    vector<type_cluster> outputPairs;
    type_cluster pair;
    dataTypeNSeq index=0;
    float normSim;
	
    while (index < nRead) //Stop reading vectors SimArray_i[1,nRead] when all reads are processed
    {
		maxSim=0;
		normSim=0;
        outputPairs.clear();
		
        for(dataTypeNSeq j=0;j<nRef;j++)
        {
            pair.sim= SimArray_[index][j];
            
            if (pair.sim>maxSim) //Update the maximum similarity
                maxSim=pair.sim;
				
				if (pair.sim>0) //Keep the pair (idRef,sim) only if the similarity is not 0
				{
					pair.idRef=j;
					outputPairs.push_back(pair);
				}
            
        }
		
		//write file only if the maximum similarity value is greater than beta
        normSim=static_cast<float> (maxSim)/norm; //normalized maximum similarity
        if (normSim>beta)
        {
            fprintf(OutFileClusters, "%.5f", normSim);
            for (std::vector<type_cluster>::iterator it = outputPairs.begin() ; it != outputPairs.end(); ++it)
            {
                normSim=static_cast<float> ((*it).sim)/norm;
                fprintf(OutFileClusters, "\t%u\t%.5f",(*it).idRef, normSim);
            }
        }
		
        //Read processed ---> new line
        fprintf(OutFileClusters,"\n");
        index++;
       
    }//end-while
	
}


int main(int argc, char **argv) {
    
    time_t t_refine=0, t_total=0;
    clock_t c_refine=0, c_total=0;

	 
	if( argc != 4)
    {
      std::cerr << "Error usage " << argv[0] << " fileFasta readLen beta" << std::endl;
       exit(1);
    }
	
    string fileFasta=argv[1];
    dataTypeSim readLen;
    float beta;
    sscanf(argv[2], "%hhu", &readLen);
    sscanf(argv[3], "%f", &beta);
	
	if ( (readLen>USim_MAX) && (dataTypeNumSim==0))
	{
		std::cerr << "Error Usage: readLen <= USim_MAX, please change settings of dataTypeNumSim in Tools.h" << std::endl;
		exit(1);
	}
	
    dataTypeNSeq numRead;
    dataTypeNSeq numRef;
    dataTypelenSeq minLCP;
    dataTypeNChar maxLen;
    string fileaux=fileFasta.substr(0,fileFasta.find(".fasta"))+".out";
    
    FILE * outAux = fopen(fileaux.c_str(), "rb");
    if (outAux==NULL) {
        std::cerr << "Error opening " << fileaux << "." << std::endl;
        printf("fopen failed, errno = %d\n", errno);
        exit (EXIT_FAILURE);
    }
	
	fread(&numRead,sizeof(dataTypeNSeq),1,outAux);
    fread(&numRef,sizeof(dataTypeNSeq),1,outAux);
    fread(&minLCP,sizeof(dataTypelenSeq),1,outAux);
    fread(&maxLen,sizeof(dataTypeNChar),1,outAux);
    
    fclose(outAux);
    
    if (maxLen>sizeMaxBuf)
        cerr << "Error Usage: maximum cluster size is " << maxLen << " greater than sizeMaxBuf, please increase sizeMaxBuf in Tools.h" << endl;
    
    dataTypelenSeq norm=readLen+1-minLCP; //to normalize similarity values
	
	//Open files .clrs, .ebwt, .da
	
	std::string fnBWT, fnCluster, fnDA;
    fnBWT = fileFasta+ ".ebwt";
    fnDA=fileFasta+".da";
    std::stringstream ssin;
    ssin << fileFasta << "." << minLCP << ".clrs\0";
	fnCluster=ssin.str();
    
	FILE * InCluster;
    InCluster = fopen(fnCluster.c_str(), "rb");
    if (InCluster==NULL) {
        std::cerr << "ClusterBWT: Error opening " << fnCluster << "." << std::endl;
        printf("fopen failed, errno = %d\n", errno);
        exit (EXIT_FAILURE);
    }
    
    fseek(InCluster, 0, SEEK_SET);


	FILE * InBWT;
    
    InBWT = fopen(fnBWT.c_str(), "rb");
    if (InBWT==NULL) {
        std::cerr << "ClusterBWT: Error opening " << fnBWT << "." << std::endl;
        printf("fopen failed, errno = %d\n", errno);
        exit (EXIT_FAILURE);
    }
    
    FILE * InDA;
    
    InDA = fopen(fnDA.c_str(), "rb");
    if (InDA==NULL) {
        std::cerr << "ClusterBWT: Error opening " << fnDA << "." << std::endl;
        printf("fopen failed, errno = %d\n", errno);
        exit (EXIT_FAILURE);
    }
  
    //Similarity arrays
    vector < vector<dataTypeSim> > SimArray_;
    SimArray_.resize(numRead);
    for(dataTypeNSeq i=0;i<numRead;i++ )
    {
        SimArray_[i].resize(numRef);
        for (dataTypeNSeq j=0;j<numRef;j++)
            SimArray_[i][j]=0;
    }
    
    time_start(&t_total, &c_total); //start time
    
    cerr << "Computing similarity arrays SimArray_i[1,numRead]..." << endl;
    
    time_start(&t_refine, &c_refine); //start time for clusterRefine
	
    dataTypeNSeq totSeq=numRead+numRef;
    dataTypeNChar clu =clusterRefine(InCluster,InDA,InBWT,numRead, totSeq, SimArray_);

	
    fclose(InCluster);
    fclose(InBWT);
    fclose(InDA);
    
    fprintf(stderr,"TIME: %.6lf\n", time_stop(t_refine, c_refine));
	
    time_start(&t_refine, &c_refine); //start time for clusterChoose
    std::stringstream ssout;
    ssout << "Clustering_results_B" << beta << "_" << fileFasta << ".txt\0";
	string fnF=ssout.str();
	//Open Output file
	FILE *fdResultF=fopen(fnF.c_str(), "w");
    if(fdResultF==NULL) {
        std::cerr << "Error opening " << fnF << "." << std::endl;
        printf("fopen failed, errno = %d\n", errno);
        exit (EXIT_FAILURE);
    }
	
    cerr << "Writing " << fnF << endl;
    
    clusterChoose(fdResultF, norm, beta, numRead, numRef, SimArray_);
	
    fclose(fdResultF);
    
    fprintf(stderr,"TIME: %.6lf\n", time_stop(t_refine, c_refine));
	
    cout << "Cluster refinement completed. Number of analyzed clusters: " << clu << "." << endl;
    cout << "Similarity arrays computed with beta=" << beta << "." << endl;
	fprintf(stdout,"Time: %.6lf\n", time_stop(t_total, c_total));
	return 0;
}

