% putRawCSC.m
% demo on how to use Mat2NlxCSC
%
%
% example:
% putRawCSC( filenameOut, timestamps, dataSamples, ['Test header info'] );
% where timestamps and dataSamples are as returned by:
% [timestamps,dataSamples] = getRawCSCData( filename, fromInd, toInd, mode );
%
%urut/MPI/dec11
function putRawCSC( filenameOut, timestamps, dataSamples,headerLine )
blocksize=512;

AppendFile=0;
ExtractMode=1;
ModeArray=1;

NumRecs=length(timestamps);

FieldSelection(1) = 1; %timestamps
FieldSelection(2) = 0; %channel nr
FieldSelection(3) = 0; %sample freq
FieldSelection(4) = 0; % valid samples
FieldSelection(5) = 1; %samples
FieldSelection(6) = 1;%header

HeaderOut{1} = ['######## Neuralynx'];     %this is REQUIRED as header prefix
HeaderOut{2} = ['FileExport Mat2NlxCSC-urut unix-vers'];    
HeaderOut{3} = [' ' headerLine];    

SamplesOut = reshape(dataSamples, blocksize, length(dataSamples)/blocksize);

Mat2NlxCSC( filenameOut, AppendFile, ExtractMode, ModeArray, NumRecs, FieldSelection, timestamps, SamplesOut, HeaderOut' );       
