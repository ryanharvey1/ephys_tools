/*---------------------------------
* Write TT
* MEX file
*
* ADR 1998
*
* input:
*    fn = file name string
*    t = n x 1: timestamps of each spike in file
*    wv = n x 4 x 32 waveforms
* 
* output: none
*
* version L4.0
--------------------------------*/

#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <matrix.h>

void WriteHeader(FILE *fp)
{
	fprintf(fp, "%%%%BEGINHEADER\n");
	fprintf(fp, "%% TT file written out by matlab.\n");
	fprintf(fp, "%% Generated by WriteTT0 written by ADR, version L4.0.\n");
	fprintf(fp, "%%%%ENDHEADER\n");
}

int bigendianMachine(void)
{
	/* returns 0 if is a littleendian machine, else returns nonzero */
	/* key is that it looks to see if short's second byte is the low order bits */
	short fullLoByte = 0xFF;
	char *byteOrder = (char *) &fullLoByte;
	return (byteOrder[1]);	
}

void mexFunction(int nOUT, mxArray *pOUT[],
				 int nINP, const mxArray *pINP[])
{
	int errorstatus;
	int bigendianFlag = bigendianMachine();
	long postHeaderPos;
	int fnlen;
	char *fn;
	FILE *fp;
	int nSpikes;
	const int *WVdim;
	double *t;
	double *wv;
	int i,j, trode, ch;  /* counters */

	/* check number of arguments: expects 3 inputs */
	if (nINP != 3)
			mexErrMsgTxt("Call with fn,t,wv as inputs."); 
	if (nOUT != 0)
		mexErrMsgTxt("Produces no output.");

	/* read inputs 1: fn */
	fnlen = (mxGetM(pINP[0]) * mxGetN(pINP[0])) + 1;
	fn = (char *) mxCalloc(fnlen, sizeof(char)); 
	if (!fn)
		mexErrMsgTxt("Not enough heap space to hold converted string.");
	errorstatus = mxGetString(pINP[0], fn,fnlen);    
	if (errorstatus)
		mexErrMsgTxt("Could not convert string data.");

	/* read inputs 2: t */
	if (mxGetN(pINP[1]) != 1)
		mexErrMsgTxt("Input 2 should be n by 1.\n");
	nSpikes = mxGetM(pINP[1]);
	t = mxGetPr(pINP[1]);

	/* read inputs 3: wv */
	if (mxGetNumberOfDimensions(pINP[2]) != 3)
		mexErrMsgTxt("Input 3 should be 3d.\n");
	WVdim = mxGetDimensions(pINP[2]);
	if (WVdim[0] != nSpikes || WVdim[1] != 4 || WVdim[2] != 32)
		mexErrMsgTxt("Input 3 should be 128 x 4 x 32.\n");
	wv = mxGetPr(pINP[2]);
		
	/* open file */
	fp = fopen(fn, "wb");
	mxFree(fn);
	if (!fp)
		mexErrMsgTxt("Could not open file.");

	mexPrintf("Writing %d Spikes\n", nSpikes);

	/* write header */
	WriteHeader(fp);

	/* write t, wv */
	for (i = 0; i < nSpikes-1; i++)			/* WHY -1?  SEEMS TO WORK BUT WHY? */
	 {
		/* write timestamp */
		union{
			char c[4];
			unsigned long ul;
		} tmpT;

		tmpT.ul = t[i];
		if (bigendianFlag)
		 {
			fwrite(&tmpT.c[0], sizeof(char), 1, fp);
			fwrite(&tmpT.c[1], sizeof(char), 1, fp);
			fwrite(&tmpT.c[2], sizeof(char), 1, fp);
			fwrite(&tmpT.c[3], sizeof(char), 1, fp);
		 }
		else /* littleendian */
		 {
			fwrite(&tmpT.c[3], sizeof(char), 1, fp);
			fwrite(&tmpT.c[2], sizeof(char), 1, fp);
			fwrite(&tmpT.c[1], sizeof(char), 1, fp);
			fwrite(&tmpT.c[0], sizeof(char), 1, fp);
		 };

		/* write wv */
		for (ch = 0; ch < 32; ch++)
		 for (trode = 0; trode < 4; trode++)
		  {
			union {
				char c[2];
				short s;
			} tmpWV;
			int subs[] = {i, trode, ch};
			int index = mxCalcSingleSubscript(pINP[2], 3, subs);  
			tmpWV.s = wv[index];
			if (bigendianFlag)
			 {
				fwrite(&tmpWV.c[0], sizeof(char), 1, fp);
				fwrite(&tmpWV.c[1], sizeof(char), 1, fp);
			 }
			else
			 {
				fwrite(&tmpWV.c[1], sizeof(char), 1, fp);
				fwrite(&tmpWV.c[0], sizeof(char), 1, fp);
			 };

		}
	 }

	fclose(fp);
}