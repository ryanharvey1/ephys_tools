/*-----------------------------------
* binning by interval
* MEX file
* 
* cowen 2002
* 
* input: 1: Spike or other timestamped data -- n_spikes x 1, assumed to be sorted
* input: 2: n_intervals x 2 matrix of intervals. The spike times will be binned into 
*			 these intervals. Also assumed to be sorted. The spikes are counted from 
the values in col 1 up to but NOT including the values in col 2.
* output: A vector of length n_intervals that contains the spike counts to each element in the bin.

  * NOTE: ASSUMES TIMESTAMPS ARE SORTED AND THE INTERVALS ARE SORTED BY THE FIRST COLUMN.
*
-----------------------------------*/

#include "mex.h"
#include <math.h>
#include <matrix.h>
//___________________________________________________________________________________
int binsearch(int n, double* data, double key, int start)
{
	// n = number of records in data.
	// data is an array of records to search
	// key is the value to find.
	// 
	// binsearch returns the index of the closest item in data
	//
	// from ADR adapted by cowen. 
	// fixed to return always the lower integreal timestamp (floor(*data)) if key is not integral (PL)
	
	int end = 0; 
	int mid = 0;
	double tmp;
	
	// binary search 
	end = n;   //-1;
	while (start < (end-1))
	{
		mid = (int) floor((double)(start + end)/2);
		tmp = floor(data[mid]);
		if (key == tmp)
			start = end = mid;
		if (key < tmp)
			end = mid;
		if (key > tmp)
			start = mid;
	}
	
	return start;
}

void mexFunction(
				 int nOUT, mxArray *pOUT[],
				 int nINP, const mxArray *pINP[])
{
	int n_intervals, n_timestamps;
	double *timestamps;
	double *intervals;
	double *result;
	int i, j,idx, last_idx, temp_time;
	/*********************************************************/
	/* check number of arguments: expects 2 inputs, 1 output */
	/*********************************************************/
	if (nINP != 2)
		mexErrMsgTxt("Call with timestamp vector and an interval matrix as inputs.");
	if (nOUT != 1)
		mexErrMsgTxt("Requires one output.");
	
	/*******************/
	/*  unpack inputs  */
	/*******************/
	n_timestamps	= mxGetNumberOfElements(pINP[0]);
	n_intervals   	= mxGetM(pINP[1]);
	
	timestamps	= (double *) mxGetPr(pINP[0]);
	intervals = (double *) mxGetPr(pINP[1]);
	
	/****************/
	/* pack outputs */
	/****************/
	pOUT[0] = mxCreateDoubleMatrix(n_intervals, 1, mxREAL);
	
	if (!pOUT[0]) 
		mexErrMsgTxt("Cannot create output array. Probably out of memory.");	
	
	result = (double *) mxGetPr(pOUT[0]);
	
	// If nothing is passed in, just return a vector of 0s.
	if (mxIsEmpty(pINP[0]) || mxIsEmpty(pINP[1]))
		return;
	
	/*********************************************************/
	/*  Do you iterate through the intervals and search for spikes for each 
	*  interval or do you iterate through the spikes and search for the proper 
	*  interval for each spike? */
	/*********************************************************/
	
	/*********************************************************/
	/*  How about this: go through each interval, find the first 
	spike within that interval, then keep going through
	the spikes until you hit one greater than the interval, 
	then move to the next interval */
	/*********************************************************/
	
	last_idx = 0;
	// Choose the most efficient method of binning depending on the data.
	//	if(n_intervals < n_timestamps)
	//if(n_intervals < n_timestamps)
	if(n_intervals < n_timestamps)
	{
	//	mexPrintf("By interval.\n");
		for (i=0;i<n_intervals;i++)
		{
			/*********************************************************/
			/* find the timestamp closest to the spike time				 */
			/*********************************************************/
			idx = binsearch(n_timestamps, timestamps, intervals[i], last_idx);
			/*********************************************************/
			/* if the timestamp is less than the start time, go to the next one */
			/*********************************************************/
			
			if (timestamps[idx] < intervals[i] && idx < n_timestamps-1)
				idx++;
			last_idx = idx;
			
			
			while(timestamps[idx] >= intervals[i] && timestamps[idx] < intervals[n_intervals + i] && idx < n_timestamps)
			{
				result[i]++;
				idx++;
			}
		}
	}else
	{
	//	mexPrintf("By timestamps.\n");
		for (i=0;i<n_timestamps;i++)
		{
			/*********************************************************/
			/* find the interval closest to the spike time				 */
			/*********************************************************/
			idx = binsearch(n_intervals, intervals, timestamps[i], last_idx);
			/*********************************************************/
			/* Move to the first interval that is near the current timestamp. */
			/*********************************************************/
			//mexPrintf("%i)  %f %f  %i ts %f \n",i,intervals[idx],intervals[idx+n_intervals],idx,timestamps[i] );
			
			while (timestamps[i] >= intervals[idx] &&  timestamps[i] < intervals[n_intervals + idx] && idx > 0)
				idx--;
			last_idx = idx;
			//idx++;
			if (idx==1)
				idx = 0;


			//mexPrintf("%i)  %f %f  %i \n",i,intervals[idx],intervals[idx+n_intervals],idx);
			
			while(timestamps[i] >= intervals[idx] && idx < n_intervals )
			{
				if(timestamps[i] >= intervals[idx] && timestamps[i] < intervals[n_intervals + idx] )
				{
					//mexPrintf("  ts %f idx %i\n",timestamps[i],idx);

					result[idx]++;
				}
				idx++;

			}
		}


	}
}
