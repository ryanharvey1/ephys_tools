% GetNlxFileProperties.dll
% 
% This function reads a file and returns some of its properties.  The function has the following format:
% 
% [NumRecs, BegTimestamp, EndTimestamp] = GetNlxFileProperties( Filename );
% 
% NumRecs - This is the number of records in the file.
% GetTimestamp - This is the first timestamp in the file.
% EndTimestamp - This is the last timestamp in the file.
% Filename - This is the filename for which we are requesting properties from.
%            This accepts any valid Neuralynx data file.
% 
% 
