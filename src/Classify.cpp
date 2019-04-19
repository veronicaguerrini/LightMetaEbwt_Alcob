#include "Tools.h"

/*
Classify ASSIGNS EACH READ TO A UNIQUE GENOME BY SIMILARITY VALUES

	This task can take in input several .txt files each one obtained from running ClusterBWT, as to consider both the forward and the reverse complement strand of single-end or paired-end reads.
	The number of .txt files (1, 2 or 4) given in input to Classify MUST be specified.

Input: The number of .txt files, list of names of .txt files obtained from ClusterBWT, total number of genomes
		n.b. all the files have the same number of lines (equal to the total number of reads) 

		The taxonomy information are provided in Datasets/Reference_database.csv
		
Output: A .txt file containing classification assignments for each read. The assignment output has 4 columns.
		The first column specifies if the read is not classified (U), classified ambiguously (A), or classified (C).
		The second column is the read IDSeq.
		The third column is either NA (in case of U or A) or the taxID/genome to which the read is assigned  
		The fourth column is the similarity score for the classification
*/

void readCorrespondence(dataTypeNSeq &numTarg, std::vector<dataTypeSet> &v_corRef) 
{
	string fileTaxID="Reference_database.csv";
	std::ifstream f_taxid;
	f_taxid.open(fileTaxID.c_str(), std::ifstream::in);
	
	if (f_taxid.is_open())   
	{
		std::cerr << "Reading taxonomical information from " << fileTaxID << ".\n";
	}
	else 
		std::cerr << "Error opening file " << fileTaxID << ".\n";
			
	string str1,str2,str3;
	getline(f_taxid,str1,';');
	getline(f_taxid,str2,';');
	getline(f_taxid,str3,'\n');
		
	getline(f_taxid,str1,';');
	getline(f_taxid,str2,';');
	getline(f_taxid,str3,'\n');
			
	while (f_taxid.good()) 
	{
		#if TaxLevel == 1
			v_corRef.push_back(atoi(str3.c_str()));
		#else
			v_corRef.push_back(str1);
		#endif
		
		getline(f_taxid,str1,';');
		getline(f_taxid,str2,';');
		getline(f_taxid,str3,'\n');
	}
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
	val=static_cast<bool> (is >> maxSim);
    while( val == 1 ) {
		val = static_cast<bool> (is >> ReadIdRef); //genome identifier index
		if ( val == 1 ) {
			val = static_cast<bool> (is >> Sim);			
			if (maxSim - Sim < (static_cast<float>(ERROR)) )
				set.insert(v_corRef[ReadIdRef]);//store only genome identifiers for which the similarity value is close to the maximum value maxSim
		}
    }
}

bool compareFirst (std::pair<float,uint> a,std::pair<float,uint> b) { return isgreater(b.first,a.first); }

void AmbiguousClassification(std::vector<string> &lineInput, uint &numFile, dataTypeNSeq idSeqRead,std::vector<dataTypeSet> &v_corRef, ofstream &out, dataTypeNSeq &numC, dataTypeNSeq &numA){
	
    std::vector<occ> CompareMates;
	occ el;
	el.TaxID=0,	el.maxRead=0, el.maxMate=0;
	
	float maxSim, Sim;
	dataTypeNSeq ReadIdRef;
				
	//Mate 1
    dataTypeNSeq index=0;
	for(uint i=0;i<numFile/2;i++ )
	{
		istringstream is(lineInput[i]);
		is >> maxSim;
		while ((is >> ReadIdRef) && (is >> Sim))
		{
			if (maxSim - Sim < (static_cast<float>(ERROR)) )
			{
				el.TaxID=v_corRef[ReadIdRef];
				el.maxRead=Sim;
				el.maxMate=0;
				index=0;
				while ((index<CompareMates.size()) && (el.TaxID <CompareMates[index].TaxID))
					index++;
						
			//Possibly update maxRead
				if ( ( index < CompareMates.size() ) && (el.TaxID == CompareMates[index].TaxID) )
                {
					if (el.maxRead > CompareMates[index].maxRead)
						CompareMates[index].maxRead=el.maxRead;
				}
				else//TaxID not in the list --> insert it
				{
					if ( index < CompareMates.size()) //index is not pointing the end of CompareMates
						CompareMates.insert(CompareMates.begin()+index,el);
					else    //append the new entry
						CompareMates.push_back(el);
                }
								
            }
        }//end-while
        
    }//end-for
	
    //If CompareMates is empty --> read classified ambiguous
    if (CompareMates.size() == 0) {
        out << "A," << idSeqRead << ","	<< "NA,0\n";
        numA++;
	return;
    }
					
    //Mate 2
	for(uint i=numFile/2;i<numFile;i++ )
	{
		istringstream is(lineInput[i]);
		is >> maxSim;
		while ((is >> ReadIdRef) && (is >> Sim))
        {
			if (maxSim - Sim < (static_cast<float>(ERROR)) )
			{
				el.TaxID=v_corRef[ReadIdRef];
				el.maxRead=0;
				el.maxMate=Sim;
                index=0;
				while ((index<CompareMates.size()) && (el.TaxID <CompareMates[index].TaxID))
                    index++;
							
                //Possibly update maxMate
                if ( ( index < CompareMates.size() ) && (el.TaxID == CompareMates[index].TaxID) ) {
                    if (el.maxMate > CompareMates[index].maxMate)
						CompareMates[index].maxMate=el.maxMate;
                }
            }
        }//end-while
    }//end-for
				 
    dataTypeNSeq TaxIDMaxSum=0;
    float maxSum=0;
    float secondMaxSum=0;
			
	for(uint i=0; i<CompareMates.size(); i++)
	{
        if (CompareMates[i].maxMate>0)
        {
            if (maxSum <= CompareMates[i].maxRead + CompareMates[i].maxMate) 
            {
                secondMaxSum = maxSum;
                maxSum = CompareMates[i].maxRead + CompareMates[i].maxMate;
                TaxIDMaxSum = CompareMates[i].TaxID;
							
            }
			else
            {
				if (secondMaxSum < CompareMates[i].maxRead + CompareMates[i].maxMate)
					secondMaxSum = CompareMates[i].maxRead + CompareMates[i].maxMate;
            }
        }
        
    }//end-for
				
	if ( maxSum > secondMaxSum + 0.001)
    {
		out << "C," << idSeqRead << "," << TaxIDMaxSum << "," << maxSum << ",0\n";
        numC++;
	}
	else
	{
		out << "A," << idSeqRead << ","	<< "NA,0\n";
        numA++;
    }
										
    CompareMates.clear();
}


	
int main(int argc, char **argv) {

	time_t t_total=0;
    clock_t c_total=0;
	
	time_start(&t_total, &c_total);

	if (argc < 2)
	{
      std::cerr << "Error usage " << argv[0] << " N fileInput1 fileInput2 ... fileInputN numGenomes fileOutput" << std::endl;
       exit(1);
    }
	
	uint numFile;
	sscanf(argv[1], "%u", &numFile);

	if( argc != numFile+4 )
    {
      std::cerr << "Error usage " << argv[0] << " N fileInput1 fileInput2 ... fileInputN numGenomes fileOutput" << std::endl;
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
	
	
	string fileOutput;
	fileOutput = argv[numFile+3];
	 
	
	std::vector<dataTypeSet> v_corRef;	//used to make correspond each reference genome identifier to a particular taxon
	cerr << "readCorrespondence: genome --> TaxID" << endl;
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
	
	std::ofstream out;
	out.open(fileOutput.c_str(), std::ofstream::out);
	if (!out.is_open())
        cerr << "ERROR: File Output not Open" << endl;
	
	out << "C/U/A,IdSeqRead,TaxID,maxSim\n";

	//Read all .txt files simultaneously
	std::set<dataTypeSet> setTarg;
	
    std::vector< std::pair<float,uint> > vector_multiset;  //vector storing in which file the max similarity value is
	std::pair<float,uint> highest_pair;
    float highest;
    uint sizeMax_multiset=0;
	bool notClass=true;
    
    cerr << "Start comparing..." << endl;
	
	for(uint i=0;i<numFile;i++ )
		getline(*fdResult_[i], lineInput[i]);
		
	
    while (!(fdResult_[0]->eof()))
	{
		for(uint i=0;i<numFile;i++ )
		{
			ReadMax(lineInput[i],&maxSim_[i]);
			if (maxSim_[i]>0)
			{
				vector_multiset.push_back(std::make_pair(maxSim_[i],i));
				notClass=false;
			}
		}
		
		//Case all the lines are empty --> read NOT CLASSIFIED
		if (notClass)
		{
			out << "U," << idSeqRead << ",NA,0\n";
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
                if(highest-curr_pair.first < static_cast<float>(ERROR) )
                    sizeMax_multiset++;
            }
            
            assert(sizeMax_multiset>0 && sizeMax_multiset<=numFile);
			
			for (uint i=0; i<sizeMax_multiset;i++)
            {
                highest_pair=*(vector_multiset.rbegin()+i);
				ReadSimilarity(lineInput[highest_pair.second], setTarg, v_corRef);
			}

			assert(setTarg.size()>0);
	
			//Assign the read to a UNIQUE TAXID --> CLASSIFIED
	
			if (setTarg.size()==1)
			{
				std::set<dataTypeSet>::iterator it = setTarg.begin();
				out << "C," << idSeqRead << "," << *it << "," << highest << "\n";
				numC++;
			}
			else	//Assign the read to MORE THAN ONE TAXID
			{
				if (numFile>1)
                    AmbiguousClassification(lineInput, numFile, idSeqRead, v_corRef, out, numC, numA);
                else
                {
                    out << "A," << idSeqRead << "," << "NA,0\n";
                    numA++;
                }
			} //end-else  assign more than one

        }//end if-else classify
        
		idSeqRead++;
		notClass=true;
		sizeMax_multiset=0;
        vector_multiset.clear();
        setTarg.clear();
		for(uint i=0;i<numFile;i++ )
			maxSim_[i]=0;
		
		for(uint i=0;i<numFile;i++ )
			getline(*fdResult_[i], lineInput[i]);
		
	}//end-while (end reading Clustering_results)

	assert(numC+numNC+numA==idSeqRead);
	for(uint i=0;i<numFile;i++ )
	{
		assert(fdResult_[i]->eof());
		fdResult_[i]->close();
	}
	
    out.close();
	
	std::cout << "Classification process completed.\nNumber of successfully classified reads: " << numC << "/" << idSeqRead << "," << std::endl;
	std::cout << "\tambiguously classified reads: " << numA << "/" << idSeqRead << "," << std::endl;
	std::cout << "\tnot classified reads: " << numNC << "/" << idSeqRead << "." << endl;
    fprintf(stdout,"Time: %.6lf\n", time_stop(t_total, c_total));
	return 0;
}
