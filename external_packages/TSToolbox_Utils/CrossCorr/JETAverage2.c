 /* ETAverage event-triggered average
  * MEX file 
  * fpbattaglia 2007 
  * input  ev: an event time series 
  *        t1 time series  in 1/10000 sec
  *               (assumed to be sorted)
  *        d the values in the dataset 
  *        binsize: the size of the bin in ms
  *        nbins: the number of bins to compute 
  * output C the event triggered average 
 *         B (optional) a vector with the times corresponding to the bins
 */







 /* adepted from PECorr original header
 * MEX file
 * 
 * batta 2002
 * 
 * input: ev: an event time series 
 *        t1, t2: two time series  in 1/10000 sec
 *                (assumed to be sorted) the correlation coefficient
 *        as a function of  time gap from eventr times will be computed
 *        binsize: the binsize for the cross corr histogram in msec
 *        nbins: the number of bins
 * output: C the peri-event correlation coefficient 
 *         B (optional) a vector with the times corresponding to the bins
 *         
 * version 0.1
 -----------------------------------*/

#define EV1_IDX 0
#define EV2_IDX 1
#define T_IDX 2
#define D_IDX 3 
#define BINSIZE_IDX 4
#define NBINS_IDX 5

#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <matrix.h>

void mexFunction(
  int nOUT, mxArray *pOUT[],
  int nINP, const mxArray *pINP[])
{
  double *ev1;
  double *ev2;
  
  double *t, *D;
  double binsize;
  double *C, *B;

  
  double dNow, w, lbound, rbound1, rbound2;
  
  int nbins, nt, nev1, nev2, nNow, *nSamples;
  int i = 0, ie2 = 0, ie1 = 0, le2, l,j,k,valk;
  int isempty = 0;


  /* check number of arguments: expects 4 inputs, 1 or 2 outputs */
  if (nINP != 6)
    mexErrMsgTxt("Call with ev1, ev2, t, d, binsize and nbins  as inputs.");
  if (nOUT != 1 && nOUT != 2)
    mexErrMsgTxt("Requires one or two outputs.");

  /* check validity of inputs */
  if (mxGetM(pINP[EV1_IDX]) != 1 && mxGetN(pINP[EV1_IDX]) != 1)
    mexErrMsgTxt("ev1 must be a row or column vector");
  if (mxGetM(pINP[EV2_IDX]) != 1 && mxGetN(pINP[EV2_IDX]) != 1)
    mexErrMsgTxt("ev2 must be a row or column vector");
  if (mxGetM(pINP[BINSIZE_IDX]) * mxGetN(pINP[BINSIZE_IDX]) != 1)
    mexErrMsgTxt("binsize must be scalar");
  if (mxGetM(pINP[NBINS_IDX]) * mxGetN(pINP[NBINS_IDX]) != 1)
    mexErrMsgTxt("nbins must be scalar");


  if (mxGetM(pINP[T_IDX]) * mxGetN(pINP[T_IDX]) == 0)
    isempty = 1;

  if(!isempty) 
    {
      if (mxGetM(pINP[T_IDX]) != 1 && mxGetN(pINP[T_IDX]) != 1)
	mexErrMsgTxt("t must be a row or column vector");
 
      if (mxGetM(pINP[D_IDX]) != 1 && mxGetN(pINP[D_IDX]) != 1)
	mexErrMsgTxt("d must be a row or column vector");

      if(mxGetM(pINP[T_IDX]) != mxGetM(pINP[D_IDX]) ||
	 mxGetN(pINP[T_IDX]) != mxGetN(pINP[D_IDX]))
	mexErrMsgTxt("t and d must have the same dims");
      
    }
  
  /* unpack inputs */
  if(!isempty)
    {
      nev1 = mxGetM(pINP[EV1_IDX]) * mxGetN(pINP[EV1_IDX]);
      ev1 = mxGetPr(pINP[EV1_IDX]);  
      nev2 = mxGetM(pINP[EV2_IDX]) * mxGetN(pINP[EV2_IDX]);
      ev2 = mxGetPr(pINP[EV2_IDX]);  

      nt = mxGetM(pINP[T_IDX]) * mxGetN(pINP[T_IDX]);
      t = mxGetPr(pINP[T_IDX]);
      D = mxGetPr(pINP[D_IDX]);

    }
	printf("nb events1 : %d; nb events2 : %d; \n",nev1,nev2);
	printf("nb time stamps : %d; first value : %f\n",nt,D[0]);
	

  binsize = mxGetScalar(pINP[BINSIZE_IDX]);
  nbins = (int)mxGetScalar(pINP[NBINS_IDX]);

  /* we want nbins to be odd */
  if ((nbins / 2) * 2 == nbins)
    nbins++;

  pOUT[0] = mxCreateDoubleMatrix(nbins, nbins, mxREAL);
  if(nOUT >= 2)
    {
      double m;
      
      pOUT[1] = mxCreateDoubleMatrix(nbins, 1, mxREAL);
      B =  mxGetPr(pOUT[1]);
      m = - binsize * ((double)nbins / 2);
      
      for(j = 0; j < nbins; j++)
	B[j] = m + j * binsize;

    }


  if(!isempty)
    {
      C = mxGetPr(pOUT[0]);
      nSamples = (int *)mxCalloc(nbins*nbins, sizeof(int));
      
      binsize *= 10;  
      /* cross correlations */
  
      w = ((double)nbins / 2) * binsize;

      for(ie1 = 0; ie1<nev1; ie1++)
	{
/*	  printf("%d\n",ie);*/
	  lbound = ev1[ie1] - w;

	  while(ev2[ie2] < lbound && ie2 < nev2 -1)
	    ie2++;
	  while(ev2[ie2-1] > lbound && ie2 > 1)
	    ie2--;

	  while(t[i] < lbound && i < nt -1)
	    i++;
	  while(t[i-1] > lbound && i > 1)
	    i--; 

	  rbound1 = lbound;
	  le2 = ie2;

	  l = i;

	  for(j = 0; j < nbins; j++)
	     {

		nNow = 0;
		dNow = 0.;

		rbound1 += binsize;

		while(t[l] < rbound1 && l < nt-1)
		{
			nNow++;
			dNow += D[l];
			l++;
	
		}

		rbound2 = lbound;
		le2 = ie2;
	
		for(k = 0; k < nbins; k++)
			{
				valk = 0;
				rbound2 += binsize;
		
				while(ev2[le2] < rbound2 && le2 < nev2-1)
				{
					/*valk++;*/
					valk = 1;
					le2++;
				}
			
				if(nNow > 0 && valk > 0)
				{
					nSamples[k+j*nbins]++;
					C[k+j*nbins] += valk*dNow/nNow;
					
				}
			}
				/*printf("\n",j);*/
	
		}

	     }

      for(j = 0; j < nbins; j++)
	for(k = 0; k < nbins; k++)
	  {
		/*printf("%f",C[k+j*nbins]);*/
		if(nSamples[k+j*nbins]>0)
			C[k+j*nbins] /= nSamples[k+j*nbins];
	  }

      mxFree((void *)nSamples);

    }
 
}
