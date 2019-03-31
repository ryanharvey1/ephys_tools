/*---------------------------------
* LoadST0_nt
* MEX file
*
* input:
*    fn = file name string
* 
* output:
*    [t, wv]
*    t = n x 1: timestamps of each spike in file
*    wv = n x 4 x 32 waveforms
*
* version 5.1
*
* Reads NT-ST (stereo-trode) files and returns a standard wv array where the 
* the ST data are in channel 1 and 2 and channels 3-4 are padded with zeros 
*(for MClust compatibility)
*
* Checks for standard Neuralynx header if present (in Cheeath versions >= 1.3)
* and automatically skips header.
*
* PL Sept 2000 
*--------------------------------*/

#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <matrix.h>

#ifdef __GNUC__
#define __int64 long long 
#endif


const int NHEADERBYTES = 16384;	


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
	
	fread(headerflag, sizeof(char), 8, fp);  
	//mexPrintf("Headerflag =  %8s \n", headerflag);
    fseek(fp,0,0);       // reset filepointer to the beginning of file
	if (strncmp(headerflag,"########",8) == 0){
		fseek(fp,NHEADERBYTES,0);  // set filepointer after byte NHEADERBYTES
	    //mexPrintf("NT-Header skipped (%d bytes)\n",NHEADERBYTES);
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


//___________________________________________________________________________________
int GetNumberOfSpikes(char* fn){
	
	// open file 
	FILE* fp = fopen(fn, "rb");
	if (!fp)
		mexErrMsgTxt("Could not open file.");

	//skip header and determine file record size
	int new_NT_format = SkipHeader(fp);
	int recSize = 260;
	if (new_NT_format) recSize = 176; 
	if (new_NT_format) SkipCheetahNTHeader(fp);
        
	// get filesize
	int postHeaderPos = ftell(fp);     // beginnig of file after header (if any)
	fseek(fp,0,2);                     // goto end of file
	int nbytes = ftell(fp) - postHeaderPos;

	int nSpikes = nbytes/recSize - 1;    /* skip last record since it may be incomplete */
	if (new_NT_format) nSpikes = nbytes/recSize; // no need to skip last record for NT_cheetah files
	mexPrintf("Reading file %s:\nRecordSize = %d,  %d spikes, %d bytes.\n",
		fn, recSize, nSpikes, nbytes);

	// cleanup
	fclose(fp);	

	return nSpikes;
}


//___________________________________________________________________________________
void ReadST(char* fn, int nSpikes, double *t, double *wv){
	
#ifdef __GNUC__
	const long long  TIMESTAMP_MAX = 0x00FFFFFFFFFFFFFFLL;  //( = 2^52-1; the largest integer fitting
#else
	const __int64 TIMESTAMP_MAX = 0x00FFFFFFFFFFFFFF;  //( = 2^52-1; the largest integer fitting
#endif
	//   in a double IEEE mantissa (=7 bytes)
	int bigendianFlag = bigendianMachine();
	
	// NT SE (SU) record
	__int64 qwTimeStamp, qwTimeStamp0;
	long dwParams[10];
	short snData[64];
	
	
	// open file 
	FILE *fp = fopen(fn, "rb");
	if (!fp)
		mexErrMsgTxt("ERROR: Could not open file.");
	SkipCheetahNTHeader(fp);     // Skip standard Neuralynx header if present (Cheetah versions >= 1.3)
	mexPrintf("NT-Header skipped (%d bytes)\n",NHEADERBYTES);
	
	
	
	// read records and convert all to double
	for (int i = 0; i < nSpikes; i++){
		
		fread(&qwTimeStamp0,  sizeof(char),   8, fp);
		fread(dwParams,      sizeof(char),  40, fp);
		fread(snData  ,      sizeof(char),  128, fp);
		
		if(bigendianFlag){
			// convert from NT(little endian) to Sun (big endian)
			qwTimeStamp = swapbytes(qwTimeStamp0);
			if(qwTimeStamp > TIMESTAMP_MAX){
				mexPrintf(" ERROR: timestamp %d in file %s is too large to fit in a double!\n",i,fn);
				mexPrintf(" Converted timestamps MAY or MAY NOT be valid - proceed with care! \n");
			}
			t[i] = (double) qwTimeStamp/100.0;
			for (int j = 0; j<2; j++)
				for(int k = 0; k<32; k++)
					wv[i + nSpikes*j + nSpikes*4*k] = (double) swapbytes(snData[j+2*k]);  
				
		} else {
			// don't convert, just copy
			qwTimeStamp = qwTimeStamp0;
			if(qwTimeStamp > TIMESTAMP_MAX){
				mexPrintf(" ERROR: timestamp %d in file %s is too large to fit in a double!\n",i,fn);
				mexPrintf(" Converted timestamps MAY or MAY NOT be valid - proceed with care! \n");
			}
			t[i] = (double) qwTimeStamp/100.0;
			for (int j = 0; j<2; j++)
				for(int k = 0; k<32; k++)
					wv[i + nSpikes*j + nSpikes*4*k] = (double) snData[j + 2*k];  
		}
		
	}
		 
   fclose(fp);	
}	



//___________________________________________________________________________________
void mexFunction(int nOUT, mxArray *pOUT[],
				 int nINP, const mxArray *pINP[])
{
	
	/* check number of arguments: expects 1 input */
	if (nINP != 1)
				mexErrMsgTxt("Call with fn as inputs.");
	if (nOUT != 2)
				mexErrMsgTxt("Requires two outputs (t, wv).");

	/* read inputs */
	int fnlen = (mxGetM(pINP[0]) * mxGetN(pINP[0])) + 1;
	char *fn = (char *) mxCalloc(fnlen, sizeof(char)); 
	if (!fn)
		mexErrMsgTxt("Not enough heap space to hold converted string.");
	int errorstatus = mxGetString(pINP[0], fn,fnlen);    
	if (errorstatus)
		mexErrMsgTxt("Could not convert string data.");
	
	int nSpikes = GetNumberOfSpikes(fn);


	// create outputs 
	pOUT[0] = mxCreateDoubleMatrix(nSpikes, 1, mxREAL);
	double *t = mxGetPr(pOUT[0]);
	
	int wvDims[] = {nSpikes, 4, 32};    // wv = nSpikes x 4 x 32   FORTRAN (column major) array
	pOUT[1] = mxCreateNumericArray(3, wvDims, mxDOUBLE_CLASS, mxREAL);
	double *wv = mxGetPr(pOUT[1]);

	// load tt file fn into t and wv arrays 
	ReadST(fn,nSpikes,t,wv);

	// cleanup
	mxFree(fn);


}

