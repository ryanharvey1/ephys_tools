/***************************************************************************
                          FindFiringField.c  -  description
                             -------------------
    begin                : Wed May 22 2002
    copyright            : (C) 2002 by Michaël Zugaro
    email                : mzugaro@andromeda.rutgers.edu
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

double xMax,yMax,threshold,*firingMap;

void FindFiringField(double x,double y,double *firingField,bool *visited)
{
  int i = (int) (x-1)*yMax+(y-1);

  if (x < 1 || x > xMax || y < 1 || y > yMax || visited[i]) return;

  visited[i] = true;
  if (firingMap[i] >= threshold)
  {
    firingField[i]=1;
    FindFiringField(x+1,y,firingField,visited);
    FindFiringField(x-1,y,firingField,visited);
    FindFiringField(x,y+1,firingField,visited);
    FindFiringField(x,y-1,firingField,visited);
  }
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  double x,y;
  double *firingField;
  bool *visited;
  int i,j;

  if (nrhs != 6)
    mexErrMsgTxt("FindFiringField: wrong number of input arguments.");
  if (nlhs != 1)
    mexErrMsgTxt("FindFiringField: wrong number of output arguments.");

  x = *mxGetPr(prhs[0]);
  y = *mxGetPr(prhs[1]);
  xMax = *mxGetPr(prhs[2]);
  yMax = *mxGetPr(prhs[3]);
  threshold = *mxGetPr(prhs[4]);
  firingMap = mxGetPr(prhs[5]);
  plhs[0] = mxCreateDoubleMatrix(yMax,xMax,mxREAL);
  firingField = mxGetPr(plhs[0]);
  memset(firingField,xMax*yMax*sizeof(double),0);

  visited = (bool*) mxCalloc((xMax+1)*(yMax+1),sizeof(bool));
  for (i = 0;i < xMax*yMax ; ++i) visited[i] = false;

  FindFiringField(x,y,firingField,visited);

  mxFree(visited);
}
