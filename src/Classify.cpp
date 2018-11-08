#include "Tools.h"

/*
Classify ASSIGNS EACH READ TO A UNIQUE GENOME BY SIMILARITY VALUES

	This task can take in input several .txt files each one obtained from running ClusterBWT, as to consider both the forward and the reverse complement strand of single-end or paired-end reads.
	The number of .txt files (1, 2 or 4) given in input to Classify MUST be specified.

Input: The number of .txt files, list of names of .txt files obtained from ClusterBWT, total number of genomes
		n.b. all the files have the same number of lines (equal to the total number of reads) 

Output: 1) A .txt file containing classified reads, precisely each line stores the pair (R_j,R_k), where 0<=R_j<numReads and 0<=R_k<numGenomes
		2) A .txt file containing ambiguous reads, where each line stores R_j,R_k,R_j,R_s,..., where 0<=R_j<numReads and 0<=R_k,R_j,R_s,...<numGenomes
		3) A .txt file containing unassigned reads, precisely each line stores the sequence identifier R_j, with  0<=R_j<numReads, for which it is not possible to determine a genome of belonging
*/

void readCorrespondence(dataTypeNSeq &numTarg, std::vector<dataTypeSet> &v_corRef) 
{
	if (TaxLevel == 0)
    {
		for (dataTypeNSeq in=0; in<numTarg; in++)
			v_corRef.push_back(in);
    }
	else
		std::cerr << "Please provide the file containing the taxonomy information." << std::endl;
	
	assert(numTarg==v_corRef.size());

}

void ReadMax(const string & s, float *maxSim){
	istringstream is( s );
	is >> (*maxSim);
}

void ReadSimilarity(const string & s, std::set<dataTypeSet> & set, std::vector<dataTypeSet> &v_corRef){
	float maxSim, Sim;
    dataTypeNSeq ReadIdRef;
	bool val;
	istringstream is( s );
	val= (is >> maxSim);
    while( val == 1 ) {
		val = (is >> ReadIdRef); //genome identifier index
		if ( val == 1 ) {
			val =(is >> Sim);			
			if (maxSim - Sim <= ERROR+0.00005)
				set.insert(v_corRef[ReadIdRef]);//store only genome identifiers for which the similarity value is close to the maximum value maxSim
		}
    }
}

void Classification(std::set<dataTypeSet> &setTarg, dataTypeNSeq idSeqRead, ofstream &out1, ofstream &out2, dataTypeNSeq &numC, dataTypeNSeq &numA){
    
        if (setTarg.size()==1) //assign the read to a unique genome
        {
            for (std::set<dataTypeSet>::iterator it = setTarg.begin(); it != setTarg.end(); ++it)
				out1 << idSeqRead << "\t" << *it << "\n";
            numC++;
        }
        else//assign the read to more than one genome
        {
            out2 << idSeqRead << "\t";
            for (std::set<dataTypeSet>::iterator it = setTarg.begin(); it != setTarg.end(); ++it)
                out2 << *it << "\t";
            out2 << "\n";
            numA++;
        }
    
}
bool compareFirst (std::pair<float,uint> a,std::pair<float,uint> b) { return (a.first<b.first); }
	
int main(int argc, char **argv) {
	
	time_t t_total=0;
    clock_t c_total=0;
	
	time_start(&t_total, &c_total);

	if (argc < 2)
	{
      std::cerr << "Error usage " << argv[0] << " N fileInput1 fileInput2 ... fileInputN numGenomes" << std::endl;
       exit(1);
    }
	
	uint numFile;
	sscanf(argv[1], "%u", &numFile);

	if( argc != numFile+3 )
    {
      std::cerr << "Error usage " << argv[0] << " N fileInput1 fileInput2 ... fileInputN numGenomes" << std::endl;
       exit(1);
    }
	
	if (numFile>4 || numFile<=0 )
	{
		std::cerr << "Error usage " << argv[0] << ": the allowed number of input files is 1 (single strand), or 2 (both strands, or paired-end reads), or 4 (both strands of paired-end reads)" << std::endl;
	}
	
	vector<string> fileInput;
    for(uint i=2;i<numFile+2;i++ )
		fileInput.push_back(argv[i]);
	
    dataTypeNSeq numTarg; 
	sscanf(argv[numFile+2], "%u", &numTarg);
	std::vector<dataTypeSet> v_corRef;	//used to make correspond each reference genome identifier to a particular taxon
	readCorrespondence(numTarg, v_corRef);//if no taxonomy information is provided in a separate file, v_corRef behaves like the identity function
	
	//Open .txt files
	std::cout << "Reading files:";
	std::vector <ifstream*> fdResult_;
	std::vector<string> lineInput;
	std::vector<float> maxSim_;//store the maximum similarity value
	for(uint i=0;i<numFile;i++ )
	{
		lineInput.push_back("");
		maxSim_.push_back(0);
		std::cout << "\n\t" << fileInput[i];
		ifstream* f = new ifstream;
		fdResult_.push_back(f);
        fdResult_[i]->open(fileInput[i].c_str());
	}
   	std::cout << "." << std::endl;
	assert(fdResult_.size()==numFile);    
	
	//Output
	dataTypeNSeq numC=0;
	dataTypeNSeq numNC=0;
	dataTypeNSeq numA=0;
	
	dataTypeNSeq idSeqRead=0;
	
	std::ofstream out0, out1, out2;
	
	out0.open("Reads_Unclassified.txt", std::ofstream::out);
    out1.open("Reads_Classified.txt", std::ofstream::out);
    out2.open("Reads_Classified_Ambiguous.txt", std::ofstream::out);
	
	//Read all .txt files simultaneously
	std::set<dataTypeSet> setTarg;
    std::vector< std::pair<float,uint> > vector_multiset;  //vector storing in which file the max similarity value is
	std::pair<float,uint> highest_pair;
    float highest;
    uint sizeMax_multiset=0;
	bool notClass=true;
    
    cerr << "Start comparing..." << endl;
    while (getline(*fdResult_[0], lineInput[0]))
	{
		ReadMax(lineInput[0],&maxSim_[0]);
		if (maxSim_[0]>0)
		{
			vector_multiset.push_back(std::make_pair(maxSim_[0],0));
			notClass=false;
		}
		
		for(uint i=1;i<numFile;i++ )
		{
			getline(*fdResult_[i], lineInput[i]);
			assert(!(fdResult_[i]->eof()));
			ReadMax(lineInput[i],&maxSim_[i]);
			if (maxSim_[i]>0)
			{
				vector_multiset.push_back(std::make_pair(maxSim_[i],i));
				notClass=false;
			}
		}
        
		//Case all the lines are empty --> read not classified
		if (notClass)
		{
			out0 << idSeqRead << "\n";
			numNC++;
		}
        else //At least one line is not empty --> read classified
        {
            std::sort (vector_multiset.begin(), vector_multiset.end(), compareFirst);

            highest_pair=*vector_multiset.rbegin();
            highest = highest_pair.first;
       
            for (std::vector<std::pair<float,uint> >::iterator it=vector_multiset.begin(); it<vector_multiset.end(); it++)
            {
                std::pair<float,uint> curr_pair=*it;
                if(highest-curr_pair.first < ERROR)
                    sizeMax_multiset++;
            }
            
            assert(sizeMax_multiset>0 && sizeMax_multiset<=numFile);
			
			for (uint i=0; i<sizeMax_multiset;i++)
            {
                highest_pair=*(vector_multiset.rbegin()+i);
				ReadSimilarity(lineInput[highest_pair.second], setTarg, v_corRef);
			}

			assert(setTarg.size()>0);
				
			Classification(setTarg, idSeqRead, out1, out2, numC, numA);
        }//end if-else classify
        
		idSeqRead++;
		notClass=true;
		sizeMax_multiset=0;
        vector_multiset.clear();
        setTarg.clear();
		for(uint i=0;i<numFile;i++ )
			maxSim_[i]=0;
		
	}//end-while (end reading Clustering_results)

	assert(numC+numNC+numA==idSeqRead);
	for(uint i=0;i<numFile;i++ )
		fdResult_[i]->close();
    out0.close();
	out1.close();
	out2.close();
	
	std::cout << "Classification process completed.\nNumber of successfully classified reads: " << numC << "/" << idSeqRead << "," << std::endl;
	std::cout << "\tambiguously classified reads: " << numA << "/" << idSeqRead << "," << std::endl;
	std::cout << "\tnot classified reads: " << numNC << "/" << idSeqRead << "." << endl;
    fprintf(stdout,"Time: %.6lf\n", time_stop(t_total, c_total));
	return 0;
}