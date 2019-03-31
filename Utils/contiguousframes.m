function [indexout]=contiguousframes(index,nframes)
% contiguousframes: takes a binary index and locates contigous sections 
% that are at least nframes long
% 
% Input: 
%           index:   binary index 
%           nframes: min number of contigous frames you wish to have 
% 
% Output: 
%           indexout: new index containing contigous frames
% 
% Ryan E Harvey 2018

index=reshape(index,[],1);

dsig=diff([0 (abs(index')>=eps) 0]);
startIndex=find(dsig>0);endIndex=find(dsig<0)-1;
stringIndex=(endIndex-startIndex+1>=nframes);
startIndex=startIndex(stringIndex);
endIndex=endIndex(stringIndex);
indices=zeros(1,max(endIndex)+1);
indices(startIndex)=1;
indices(endIndex+1)=indices(endIndex+1)-1;
indices=find(cumsum(indices));
index=zeros(length(index),1);
index(indices',1)=1;

indexout=logical(index);
end