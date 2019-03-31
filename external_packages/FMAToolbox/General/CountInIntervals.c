/***************************************************************************
                            CountInIntervals.c
                            --------------------
    begin                : December 2009
    copyright            : (C) 2009 by Michaël Zugaro
    email                : michael.zugaro@college-de-france.fr
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "mex.h"

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  double *t,*intervals,*counts;
  int i,j,m,n;

	if (nrhs != 2) mexErrMsgTxt("Incorrect number of parameters (type 'help <a href=\"matlab:help CountInIntervals\">CountInIntervals</a>' for details).");

	intervals = mxGetPr(prhs[1]);
	if ( mxGetN(prhs[1]) != 2 ) mexErrMsgTxt("Incorrect intervals (type 'help <a href=\"matlab:help CountInIntervals\">CountInIntervals</a>' for details).");
	m = mxGetM(prhs[1]);

	t = mxGetPr(prhs[0]);
	n = mxGetM(prhs[0]);
	if ( n == 1 ) n = mxGetN(prhs[0]);

	plhs[0] = mxCreateDoubleMatrix(m,1,mxREAL);
	counts = mxGetPr(plhs[0]);

	i = j = 0;
	while ( j < n )
	{
		/* Find interval to which t[j] belongs */
		while ( t[j] > intervals[m+i] && i < m ) ++i;
		if ( t[j] > intervals[i+m] ) break;
		/* How many of t[j], t[j+1]... belong to this interval? */
		while ( t[j] <= intervals[m+i] && j < n ) { ++counts[i]; ++j; }
	}
}
