 %IMPORT OPEN FIELD DATA

%ANNOTATE  - make sure to explain what and why for a given line command
cd 'd:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\PAE_OF\'
[~, ~, raw] = xlsread('OFmax_min.xlsx');
stringVectors = string(raw(:,1));
stringVectors(ismissing(stringVectors)) = '';

raw(1,:)=[];
params=table;

for i = 1:length(raw) %"i" is an itteration and allows the loop to build
    %Get Path to Data 
    params.VideoName{i}=raw(i,1);
    params.dia(i)=202;
    params.group(i)=raw(i,6);
end

cd 'd:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\PAE_OF\TrackData\'
for i=1:length(raw)
    tempCoords = importFile([raw{i,1},'.txt'], 2, inf);
    params.rawcoords{i}=tempCoords;
end

%GET MAX/MIN
tempX=[];
tempY=[];
for ii=1:length(params.rawcoords)
    tempX=[tempX;params.rawcoords{ii}(:,1)];
    tempY=[tempY;params.rawcoords{ii}(:,2)];
end

params.Xmax(:,1)=max(tempX); params.Xmin(:,1)=min(tempX); params.Ymin(:,1)=min(tempY); params.Ymax(:,1)=max(tempY);

