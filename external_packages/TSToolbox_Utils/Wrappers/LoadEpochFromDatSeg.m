function ep = LoadEpochFromDatSeg(ix)

% ep = LoadEpochFromDatSeg(ix)

if exist('RecordingDuration.txt','file')
    datSegments = load('RecordingDuration.txt');
    datSegments = intervalSet([0;cumsum(datSegments(1:end-1))]',[cumsum(datSegments)-1]);    
elseif exist('DatSegments.txt','file')
    datSegments = load('DatSegments.txt');
    datSegments = intervalSet([0;cumsum(datSegments(1:end-1))]/1250',[cumsum(datSegments)-1]/1250);    
else 
    error('No Recoding duration file!')
end

ep = subset(datSegments,ix);