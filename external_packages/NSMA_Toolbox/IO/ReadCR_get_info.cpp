/*---------------------------------
* ReadCR_nt_partial_load
* Read a CR file from Cheetah (NT and unix Version 2.x)
* MEX file
*
*
* input:
*    fn = file name string

* output:
*    [t, cr, SampFreq]
*    t = n x 1: timestamp of each CR Record in file
*    cr = n x 512   one CR record of 512 ADC samples
*    SampFreq = 1 x 1 common Sampling Freqency of CR data (Hz) 
*
* version 2.0
* PL Sept 1999
* cowen Oct 2000 - modified for partial load
* cowen Feb 2002 - modified to load old v1.? UNIX cheetah data.
--------------------------------*/

#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <matrix.h>


#ifdef __GNUC__
#define __int64 long long 
#endif

/***********************************************************/
/* DEFINTIONS: */
/***********************************************************/


typedef struct block_type_tag {
	FILE *fp;
	double	TimeStamp;	
	long	ChannelNum;	
	long	SampleFreq;	
	long	NumValidSamples;	
	short	Data[512];		
} BLOCK_TYPE;

//___________________________________________________________________________________
int SkipCheetahNTHeader(FILE *fp)
{
	// Cheetah NT versions 1.3 or greater have a standard Neuralynx header of ascii strings of total
	// length of 16384 bytes in the beginning of the binary file. 
	// Check if Cheetah header present (in versions > 1.3) and skip header
	//
	// returns 1 if new Neuralynx header is present and was skipped
	//         0 if NO  header is present

    char headerflag[8];   // if first 8 bytes of AnyCheetahfile.dat contain "########" there is a header 
                          // of length 16,384 bytes to skip (from the very beginning including the eight '#') 
    const int NHEADERBYTES = 16384;	
	
	fread(headerflag, sizeof(char), 8, fp);  
	//mexPrintf("Headerflag =  %8s \n", headerflag);
    fseek(fp,0,0);       // reset filepointer to the beginning of file
	if (strncmp(headerflag,"########",8) == 0){
		fseek(fp,NHEADERBYTES,0);  // set filepointer after byte NHEADERBYTES
	    mexPrintf("NT-Header skipped (%d bytes)\n",NHEADERBYTES);
		return 1;
	} 
	return 0;
}

//___________________________________________________________________________________
int SkipHeader(FILE *fp)
{
	/* returns 0 if header present and skipped  (success) 
	**         1 if NO header present in file (indicates a new NT_cheetah TT file) */
	long curpos = ftell(fp);
	char headerline[81];

	fgets(headerline, 80, fp);
	if (strncmp(headerline, "%%BEGINHEADER", 13) == 0){
		while (strncmp(headerline, "%%ENDHEADER",11) != 0)
			fgets(headerline, 80, fp);
		return 0;
	} else {
		fseek(fp, curpos, SEEK_SET);
		return 1;
	}
}

//___________________________________________________________________________________
int bigendianMachine(void)
{
	/* returns 0 if is a littleendian machine, else returns nonzero */
	/* key is that it looks to see if short's second byte is the low order bits */
	short fullLoByte = 0xFF;
	char *byteOrder = (char *) &fullLoByte;
	return (byteOrder[1]);	
}


//___________________________________________________________________________________
inline short swapbytes(short ii)
// swap byte order of a short: (0,1) -> (1,0)
{
	union {
		short s;
		char c[2];
	} tmp0,tmp;
	
	tmp.s = ii;
	tmp0.c[0] = tmp.c[1];					
	tmp0.c[1] = tmp.c[0];

	return tmp0.s;
}

//___________________________________________________________________________________
inline unsigned short swapbytes(unsigned short ii)
// swap byte order of a short: (0,1) -> (1,0)
{
	union {
		unsigned short us;
		char c[2];
	} tmp0,tmp;
	
	tmp.us = ii;
	tmp0.c[0] = tmp.c[1];					
	tmp0.c[1] = tmp.c[0];

	return tmp0.us;
}


//___________________________________________________________________________________
inline unsigned long swapbytes(unsigned long ii)
// swap byte order of an unsigned long: (0,1,2,3) -> (3,2,1,0)
{
	union {
		unsigned long ul;
		char c[4];
	} tmp0,tmp;
	
	tmp.ul = ii;
	tmp0.c[0] = tmp.c[3];					
	tmp0.c[1] = tmp.c[2];
	tmp0.c[2] = tmp.c[1];					
	tmp0.c[3] = tmp.c[0];

	return tmp0.ul;
}

//___________________________________________________________________________________
inline long swapbytes(long ii)
// swap byte order of a long: (0,1,2,3) -> (3,2,1,0)
{
	union {
		long l;
		char c[4];
	} tmp0,tmp;
	
	tmp.l = ii;
	tmp0.c[0] = tmp.c[3];					
	tmp0.c[1] = tmp.c[2];
	tmp0.c[2] = tmp.c[1];					
	tmp0.c[3] = tmp.c[0];

	return tmp0.l;
}

inline long doubleswapbytes(long ii)
// swap byte order of a long: (0,1,2,3) -> (2,3,1,0)
{
	union {
		long l;
		char c[4];
	} tmp0,tmp;
	
	tmp.l = ii;
	tmp0.c[0] = tmp.c[2];					
	tmp0.c[1] = tmp.c[3];
	tmp0.c[2] = tmp.c[1];					
	tmp0.c[3] = tmp.c[0];

	return tmp0.l;
}



//___________________________________________________________________________________
inline __int64 swapbytes(__int64 ii)
// swap byte order of a long long: (0,1,2,3,4,5,6,7) -> (7,6,5,4,3,2,1,0)
{
	union {
		__int64 ll;
		char c[8];
	} tmp0,tmp;
	
	tmp.ll = ii;
	tmp0.c[0] = tmp.c[7];					
	tmp0.c[1] = tmp.c[6];
	tmp0.c[2] = tmp.c[5];					
	tmp0.c[3] = tmp.c[4];
	tmp0.c[4] = tmp.c[3];					
	tmp0.c[5] = tmp.c[2];
	tmp0.c[6] = tmp.c[1];					
	tmp0.c[7] = tmp.c[0];

	return tmp0.ll;
}



int get_data_v1(FILE *fp, BLOCK_TYPE *blk)
{	

/* For old UNIX cheetah data */

#ifdef __GNUC__
	const long long  TIMESTAMP_MAX = 0x00FFFFFFFFFFFFFFLL;  //( = 2^52-1; the largest integer fitting
#else
	const __int64 TIMESTAMP_MAX = 0x00FFFFFFFFFFFFFF;  //( = 2^52-1; the largest integer fitting
#endif
	unsigned long qwTimeStamp;
	unsigned long dwChannelNum = 999;
	unsigned long dwSampleFreq = 0;
	unsigned short dwNumValidSamples;
	int j;

	fread(&qwTimeStamp,  sizeof( long),   1, fp);
	fread(&dwNumValidSamples, sizeof( short), 1, fp);
	fread(&dwSampleFreq, sizeof( long),   1, fp);
	fread(blk->Data, sizeof( short), 512, fp);
	qwTimeStamp = swapbytes(qwTimeStamp);
	
	if(qwTimeStamp > TIMESTAMP_MAX){
		mexPrintf(" ERROR: timestamp %d is too large to fit in a double!\n",qwTimeStamp);
		mexPrintf(" Converted timestamps MAY or MAY NOT be valid - proceed with care! \n");
	}
	dwSampleFreq = swapbytes(dwSampleFreq);
	dwNumValidSamples = swapbytes(dwNumValidSamples);
	for (j = 0; j<512; j++)
		blk->Data[j] = swapbytes(blk->Data[j]);
	/* Convert into the standard block format */

	blk->TimeStamp       = (double) qwTimeStamp; /* Timestamps already in 1/10000sec */
	blk->ChannelNum      = (unsigned long) dwChannelNum ;
	blk->SampleFreq      = (unsigned long) dwSampleFreq;
	blk->NumValidSamples = (unsigned long) dwNumValidSamples;

	return (1);
}


int get_data_v2(FILE *fp, int bigendianFlag, BLOCK_TYPE *blk)
{	
	/* For newer NT cheetah data */
#ifdef __GNUC__
	const long long  TIMESTAMP_MAX = 0x00FFFFFFFFFFFFFFLL;  //( = 2^52-1; the largest integer fitting
#else
	const __int64 TIMESTAMP_MAX = 0x00FFFFFFFFFFFFFF;  //( = 2^52-1; the largest integer fitting
#endif
	__int64 qwTimeStamp;
	unsigned long dwChannelNum;
	unsigned long dwSampleFreq;
	unsigned long dwNumValidSamples;
	//signed short snSamples[512];
	int j;
	fread(&qwTimeStamp,  sizeof(__int64),   1, fp);
	fread(&dwChannelNum, sizeof(long),   1, fp);
	fread(&dwSampleFreq, sizeof(long),   1, fp);
	fread(&dwNumValidSamples, sizeof(long), 1, fp);
	fread(blk->Data, sizeof(short), 512, fp);
	/* fread(&junk,sizeof(long),1,fp); */
			
	if(bigendianFlag){
	// convert from NT(little endian) to Sun (big endian)
		qwTimeStamp = swapbytes(qwTimeStamp);
		if(qwTimeStamp > TIMESTAMP_MAX){
			mexPrintf(" ERROR: timestamp %d is too large to fit in a double!\n",qwTimeStamp);
			mexPrintf(" Converted timestamps MAY or MAY NOT be valid - proceed with care! \n");
		}
		dwChannelNum = swapbytes(dwChannelNum);
		dwSampleFreq = swapbytes(dwSampleFreq);
		dwNumValidSamples = swapbytes(dwNumValidSamples);
		for (j = 0; j<512; j++)
			blk->Data[j] = swapbytes(blk->Data[j]);
	}
	/* Convert into the standard block format */

	blk->TimeStamp = (double) qwTimeStamp/100; /* Convert timestamps to what we are used to (1/10000 sec) */
	blk->ChannelNum = dwChannelNum ;
	blk->SampleFreq = dwSampleFreq;
	blk->NumValidSamples = dwNumValidSamples;
	return (1);
}


//___________________________________________________________________________________
void mexFunction(int nOUT, mxArray *pOUT[],
					  int nINP, const mxArray *pINP[])
{
	int errorstatus;
	int bigendianFlag = bigendianMachine();
	long postHeaderPos,endPos,start_pos,end_pos;
	int fnlen;
	int new_NT_format;     /* flag for new NT_format TT files (0=old SUN, 1=new NT)  */
	char *fn;
	FILE *fp;
	int nSamples = 0;      /* to be implemented later */
	int nRecords;
	const int NEW_CR_RECSIZE = 8+sizeof(int)+2*sizeof(long)+512*sizeof(short);

	
	unsigned long junk = 0; 
	double *t;
	double *cr;
	double sampFreq, sampFreq0 = 0.0; 
	double *crptr;
	
	double start_ts = 0; // Start ts of the data to be loaded
	double end_ts = 9999999999;   // End ts of the data to be loaded
	double first_ts = 0; // First ts of the file.
	double actual_start_ts = 0; // The actual ts of the block closest to the start ts.
	double actual_end_ts = 0;
	
	BLOCK_TYPE blk; /* The data from one block in the CR file */
	int crDims[] = {nSamples, 512};
	int subs[] = {0, 0};
	int index; 
	int i,j,k,rec_count = 0;  /* counters */
	int partial_load = 0, found_start = 0, found_end = 0;     /* Flags */
	int got_data = 0;		
	double *M;
	
			
	/* check number of arguments: expects 1 input */
	if (!(nINP == 1))
		mexErrMsgTxt("Call with fn as inputs.");
	if (nOUT != 4)
		mexErrMsgTxt("Requires four outputs .");
	
	/* read inputs */
	fnlen = (mxGetM(pINP[0]) * mxGetN(pINP[0])) + 1;

	fn = (char *) mxCalloc(fnlen, sizeof(char)); 
	if (!fn)
		mexErrMsgTxt("Not enough heap space to hold converted string.");
	errorstatus = mxGetString(pINP[0], fn,fnlen);    
	if (errorstatus)
		mexErrMsgTxt("Could not convert string data.");
	
	/* open file */
	fp = fopen(fn, "rb");
	if (!fp)
		mexErrMsgTxt("Could not open file.");
	
	/* skip header */
	new_NT_format = SkipHeader(fp);
	if (new_NT_format) SkipCheetahNTHeader(fp);     // Skip standard Neuralynx header if present (Cheetah versions >= 1.3)
	postHeaderPos = ftell(fp);
	
	/* count number of Records */
	fseek(fp, 0, SEEK_END);	
	endPos = ftell(fp);
	nRecords = (int)floor((double)(endPos - postHeaderPos) / NEW_CR_RECSIZE);

	if (nSamples == 0 || nSamples > nRecords) nSamples = nRecords;
	mexPrintf("File contains %d CR records.\n", nSamples);
	
	/* read records and convert all to double*/
	fseek(fp, postHeaderPos, SEEK_SET);
	start_pos = postHeaderPos;
	end_pos = endPos;
	i = 0;
	
	while ( i < nSamples )
	{
   
		if (new_NT_format)
			got_data = get_data_v2(fp, bigendianFlag, &blk);
		else
			got_data = get_data_v1(fp, &blk);
		
		/* Find the start and end positions in the file */
		//mexPrintf(" 3Timestamp %g vs %d\n",blk.TimeStamp, blk.NumValidSamples );	

		if (!found_start && blk.TimeStamp >= start_ts)
		{
			actual_start_ts = blk.TimeStamp;
			found_start = 1;
			start_pos = ftell(fp);
		}
		if (!found_end && blk.TimeStamp>= end_ts) {
			actual_end_ts = blk.TimeStamp;
			found_end = 1;
			i = nSamples; // End
			end_pos = ftell(fp);
		}

		sampFreq = blk.SampleFreq;
		if (i>0 && sampFreq != sampFreq0){ 
			mexPrintf("New Warning: Sampling Frequency changed from %g to %g within file!!\n", 
			sampFreq0,sampFreq);   
			sampFreq0 = sampFreq;
		}
		if (i<1) {
			sampFreq0 = sampFreq;
			first_ts = blk.TimeStamp;

			mexPrintf("ts[%d]= %g \n", i,first_ts);  
			mexPrintf("ChanNum[%d]= %d \n", i,blk.ChannelNum);  
			mexPrintf("sampFreq[%d]= %d \n", i,blk.SampleFreq);  
			mexPrintf("dwNumValid[%d] = %d\n", i,blk.NumValidSamples);  
			/* mexPrintf("junk[%d] = %d\n", i,junk);  */
			mexPrintf("Sampling Frequency in %s is %g \n", 
				fn,sampFreq);  
			}
		
		// Count the number of records within the start and end range.
		if (found_start && !found_end) {
			rec_count++;
		}
		i++;
	}

	if (partial_load)
		nSamples = rec_count;

	if (actual_start_ts == actual_end_ts)
	{
		mexPrintf("Specified range not found within the data.");	
		nSamples = 0;
	}

	if (!found_end)
	{
		actual_end_ts = blk.TimeStamp;
		end_pos = ftell(fp);
	}


	// Allocate memory to the t and cr varibles.

	
	mexPrintf("Start ts: %g  End ts: %g \n", actual_start_ts, actual_end_ts);
	//mexPrintf("File Start ts: %g   File End ts: %g \n", first_ts, blk.TimeStamp);

	t  = (double *) mxCalloc(nSamples, sizeof(double));
	cr = (double *) mxCalloc(nSamples * 512, sizeof(double));
	if (!t || !cr)
		mexErrMsgTxt("Not enough heap space for TT components.");


	mxFree(cr);
	mexPrintf("Read %i records. \n", nSamples);
	/* create outputs */
	pOUT[0] = mxCreateDoubleMatrix(1,1,mxREAL);
	M = mxGetPr(pOUT[0]); M[0] = (double)actual_start_ts;
	pOUT[1] = mxCreateDoubleMatrix(1,1,mxREAL);
	M = mxGetPr(pOUT[1]); M[0] = (double)actual_end_ts;
	pOUT[2] = mxCreateDoubleMatrix(1,1,mxREAL);
	M = mxGetPr(pOUT[2]); M[0] = (double)sampFreq;
	pOUT[3] = mxCreateDoubleMatrix(1,1,mxREAL);
	M = mxGetPr(pOUT[3]); M[0] = (double)rec_count;
	
	mxFree(fn);
	
	fclose(fp);
	
}
