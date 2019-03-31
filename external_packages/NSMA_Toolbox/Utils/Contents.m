% Contents of Utils folder
% ADR 1998, version L4.0, last modified 4/03 by MRN
%
% File routines
%   FindFiles    - Finds all files in subdirectory tree that match a wildcard input globfn
%   popdir       - Pops top directory off current directory stack, cd's to it
%   pushdir      - Pushes current dir onto directory stack, cd's to new dir if given
%   ReadFileList - Reads list of files from ascii file, creates cell array of filenames
%
% Graphical routines
%   AxisEnlargeSelf          - Allows a subplot to be clicked on to be expanded
%   DisplayProgress          - Progress bar
%   DrawConvexHull           - Allows user to draw convex hull on current axis; returns x,y points on hull
%   MultiPlotCellArray       - Creates multiple plots by passing multiple sets of inputs to a given function
%   ParentFigureHandle       - Gets figure handle above current object
%   Rotate3dPlot             - Rotates current 3D figure around the z axis
%   Subplots                 - Chooses dimensions of nicely laid out grid for n subplots to be arranged in
%   TransferBetweenListboxes - Transfers values between two listboxes
%
% General utils
%   binsearch                 - (MEX FILE) Does a binary search to find index of specified value in a sorted list of data
%   Extract_varargin          - (NOT A FUNCTION) Creates variables with names and values provided in varargin list of function
%   FillWithLastNonzero       - Replaces all zero elements in a vector with the closest previous non-zero value
%   hms2ts                    - Converts time in hours,minutes,seconds into timestamp (0.1 msec) units
%   InConvexHull              - Determines whether specified points are within a specified convex hull
%   ndhist                    - (MEX FILE) Creates n-dimensional histogram
%   ReadScopeFile             - Reads propriety binary file from Nicolet 310 Storage Oscilloscope
%   SelectAlongFirstDimension - Restricts data of n-dimensional array to specified indices in the 1st dimension
%   sortcell                  - Sorts a cell array of strings
%   StatusWarning             - Prints a status warning
%   ts2hms                    - Converts time in timestamp (0.1 msec) units into hours,minutes,seconds

