% * ReadSE  
% * MEX file
% *
% * input:
% *    fn = file name string
% *    records_to_get = an array that is either a range of values
% *    record_units = 1: timestamp list.(a vector of timestamps to load (uses binsearch to find them in file))
% *                   2: record number list (a vector of records to load)
% *					3: range of timestamps (a vector with 2 elements: a start and an end timestamp)
% *					4: range of records (a vector with 2 elements: a start and an end record number)
% *
% *	 if only fn is passed in then the entire file is opened.
% *  if only fn is passed AND only t is provided as the output, then all 
% *       of the timestamps for the entire file are returned.
% *
% * output:
% *    [t, wv]
% *    t = n x 1: timestamps of each spike in file
% *    wv = n x 32 waveforms
% *
% * version 5.1
% *
% * Reads  sun-TT files and NT-SE files and distinguishes them
% * by checking if a header exists (for sun TTfiles) or not (for NT-SEfiles)
% *
% * Checks for standard Neuralynx header if present (in Cheeath versions >= 1.3)
% * and automatically skips header.
% *
% *
% * TO DO: Do the binsearch on the file and not on the entire array of timestamps.
