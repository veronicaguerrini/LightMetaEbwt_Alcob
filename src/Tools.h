/*
 * Basic tools and defines
 */

#ifndef _TOOLS_H_
#define _TOOLS_H_

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <errno.h>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <stdio.h>
#include <assert.h>
#include <sys/stat.h>
#include <math.h>
#include <numeric>
#include <malloc.h>
#include <time.h>


#define BUFFERLCPSIZE 1073741824

#define sizeMaxBuf 4096
#define BUFFERCLUSTER 1048576

#define ERROR 0.02

using namespace std;

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long ulong;

#define dataTypedimAlpha uchar  //---> eBWT (size of the alphabet) -- biologic case 6 symbols ($,A,C,G,N,T)
#define dataTypelenSeq uint     //---> LCP (length of the sequences)
#define dataTypeNumSeq 0        //---> Number of sequences in the input fasta file. (If =0 -> uint)
#define dataTypeNumChar 1       //---> Length of the list of sorted suffixes, i.e. number of characters in the input fasta file. (If =0 -> uint)
#define dataTypeNumSim 0		//---> Maximum similarity value depends on read length (If =0-> uchar)

#if dataTypeNumSeq == 1
#   define dataTypeNSeq ulong
#else
#   define dataTypeNSeq uint
#endif

#if dataTypeNumChar == 1
#   define dataTypeNChar ulong
#else
#   define dataTypeNChar uint
#endif

#if dataTypeNumSim == 1
#   define dataTypeSim uint
#else
#   define dataTypeSim uchar  
#endif

#if dataTypeNumSim == 1
#   define USim_MAX UINT_MAX
#else
#   define USim_MAX UCHAR_MAX 
#endif

struct ElementCluster {
	dataTypeNChar pStart;
	dataTypeNChar len;
};

struct ElementInCluster {
	dataTypeNSeq reads;	
	dataTypeNSeq targs; 
};

struct type_cluster {
	dataTypeNSeq idRef;		
	dataTypeSim sim;
};

#define TaxLevel 1	//Taxonomy rank

#if TaxLevel == 0
#   define dataTypeSet string
#else
#   define dataTypeSet dataTypeNSeq
#endif

struct occ {
	dataTypeSet TaxID;
	float maxRead;
	float maxMate;
};

void time_start(time_t *t_time, clock_t *c_clock){
    
	*t_time = time(NULL);
	*c_clock =  clock();
}

double time_stop(time_t t_time, clock_t c_clock){
    
	double aux1 = (clock() - c_clock) / (double)(CLOCKS_PER_SEC);
	double aux2 = difftime (time(NULL),t_time);
	
	return aux1;
}
#endif
