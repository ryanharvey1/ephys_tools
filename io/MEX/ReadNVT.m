function [Data]=ReadNVT(fn)
%[Data]=ReadNVT(fn)
%Reads in a NeuroLynx NVT file (with filename fn) into matlab
%this assumes that the format has a 16Kb ascii header and then a repeating
%binary record given by the following format
%variablename, variabletype, number_in_a_row
%RecordFormat={
% 'swstx','uint16',1;
% 'swid','uint16',1;
% 'swdata_size','uint16',1;
% 'TimeStamp','uint64',1;
% 'dwPoints','uint32',400;
% 'sncrc','int16',1;
% 'Xloc','int32',1;
% 'Yloc','int32',1;
% 'Angle','int32',1;
% 'dntargets','int32',50;
% };

%check to make sure that this file exists
fileinfo=dir(fn);
if (isempty(fileinfo))
    disp('file not found');
    output=[];
    return;
end

%open the file
fid=fopen(fn,'r','ieee-le');

%define the size of the header (in this case 16Kb)
startpoint=16*1024;

%define the format of each record
RecordFormat=define_RecordFormat;
%read in the data
Data=ReadData(fid,startpoint,RecordFormat);


function sectionsize=EvalSectionSize(fid,Format)
%EvalSectionSize(fid,Format)
%function to calculate the size of each record given the format

%convert the format to a struct
s=cell2struct(Format,{'name','numType','number'},2);

%start at the beginning
fseek(fid,0,-1);
%read in one record piece by piece
for i=1:length(s)
 eval(['[Section.' s(i).name ',n]=fread(fid,' num2str(s(i).number) ',''' s(i).numType ''');']);
end
%read out where we are in the file now, this should be the recordsize
sectionsize=ftell(fid);


function RecordFormat=define_RecordFormat
%define the format for the repeating record
%this is the format for the NVT file given by the neurolynx documentation
%should go 'variablename','variabletype',numinarow
%for single values use 1 for numinarow
%to read in arrays, use the length of the array
RecordFormat={
 'swstx','uint16',1;
 'swid','uint16',1;
 'swdata_size','uint16',1;
 'TimeStamp','uint64',1;
 'dwPoints','uint32',400;
 'sncrc','int16',1;
 'Xloc','int32',1;
 'Yloc','int32',1;
 'Angle','int32',1;
 'dntargets','int32',50;
 };

function Data=ReadData(fid,offset,Format)
%re
recordsize=EvalSectionSize(fid,Format);
s=cell2struct(Format,{'name','numType','number'},2);
cumoffset=0;
Data=[];
for i=1:length(s)
    fseek(fid,offset+cumoffset,'bof');
    eval(['dummy=' s(i).numType '(1);']);
    dummyinfo=whos('dummy');
    cumoffset=cumoffset+(dummyinfo.bytes*s(i).number);
    
 eval(['[Data.' s(i).name ']=fread(fid,''' num2str(s(i).number) '*' s(i).numType '=>' s(i).numType ''',' ...
     num2str(recordsize-(dummyinfo.bytes*s(i).number)) ');']);

  if (s(i).number>1)
      eval(['[Data.' s(i).name ']=reshape(Data.' s(i).name ',' num2str(s(i).number) ',length(Data.' s(i).name ')/' num2str(s(i).number) ')'';']);
  end
end

