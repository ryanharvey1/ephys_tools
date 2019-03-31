% ImportIgorOutput
% imports data from RH Place Cell Data Worksheet.4.xlsx to extract the original igor outputs
% packs everything neatly away in 3d matrices 
% 
% 
% 
%% LOAD WORKSHEETS
clear;close all;clc


data=importdata('/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/RH Place Cell Data Worksheet.4.xlsx');

%% FORMAT SESSION ID TO MATCH WITH PLACE CELL SESSION ID BELOW
ratIDs=fieldnames(data.data);
ratIDs=ratIDs(contains(ratIDs,'Cave'));
for i=1:length(ratIDs)
    % extract id column
    temptxt=data.textdata.(ratIDs{i});
    id=temptxt(:,1);
    id=id(3:end,1);
    
    % extract tetrode, cell number, and date
    tempdata=data.data.(ratIDs{i});
    tempdata=tempdata(:,1:3);
    
    % check to see if the TT number is there in tempdata
    if sum(~isnan(tempdata(:,2)))==0
        B = regexp(temptxt(3:end,3),'\d*','Match');
        for ii=1:length(B)
            if ~cellfun(@isnumeric, B{ii,1})
                B2(ii,1)=str2double(B{ii,1});
            end
        end
        B2(B2==0)=NaN;
        B2=[B2;NaN;NaN;NaN;NaN];
        tempdata(:,2)=B2;
    end
    
    % format date
    for ii=1:length(tempdata)
        if ~isnan(tempdata(ii,1))
        tempdate=datestr(tempdata(ii,1),2);
        tempdate=strsplit(tempdate,'/');
        date{ii,1}=[['20',tempdate{3}],' ',[tempdate{1},'-',num2str(str2double(tempdate{2})-1)],' ',id{ii}];
        else
%             tempdate{ii,1}=NaN;
        end
    end
    clear tempdate ii
    
    % recreate cell numbers if they don't exist
    if sum(~isnan(tempdata(:,3)))==0
        TT=tempdata(:,2);
        for ii=1:length(TT)
           if ii==1
               cell(ii,1)=1;
           elseif isnan(TT(ii))
               cell(ii,1)=NaN;
           elseif TT(ii)==TT(ii-5)
               cell(ii,1)=cell(ii-5)+1;
           elseif TT(ii)~=TT(ii-5)
               cell(ii,1)=1;
           end
        end
        tempdata(:,3)=cells;
    end
    
    for ii=1:length(date)
        paths{ii,1}=[date{ii},'/ASCII/TT',num2str(tempdata(ii,2)),' SS_',num2str(tempdata(ii,3))];
    end
    
    for ii=1:length(paths)
        if contains(paths{ii},'NaN')==1
            paths{ii,1}=NaN;
        end
    end
    paths=[paths;{NaN};{NaN};{NaN};{NaN}];
   
    ALL_paths{i}=paths;
    clearvars -except ratIDs i ALL_paths data
end
everypath=[ALL_paths{1};ALL_paths{2};ALL_paths{3};ALL_paths{4};ALL_paths{5};ALL_paths{6};ALL_paths{7};ALL_paths{8};ALL_paths{9};ALL_paths{10};ALL_paths{11}];
everypath(cellfun(@(x) isnumeric(x) && isnan(x), everypath)) = {''};

%% GET PATH NAMES TO PLACE CELLS
load('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/GW_Field_StatsPlaceCells_Tilted_Mice.mat')
controlid=data.textdata.Field_StatsPlaceCells0x2DTilted(ismember(data.textdata.Field_StatsPlaceCells0x2DTilted(:,2),'Control'),1);
tiltedid=data.textdata.Field_StatsPlaceCells0x2DTilted(ismember(data.textdata.Field_StatsPlaceCells0x2DTilted(:,2),'Tilted'),1);

% get path names for first session
controlid=controlid(1:5:end,:);
tiltedid=tiltedid(1:5:end,:);

controlformatedpaths=formatpaths(controlid);
tiltedformatedpaths=formatpaths(tiltedid);

%% GET INDEX OF PLACE CELLS
for i=1:length(controlid)
    Index(:,i) = contains(everypath,controlformatedpaths{i},'IgnoreCase',true);
end
controlIndex=sum(Index,2);
clear Index

for i=1:length(tiltedid)
    Index(:,i) = contains(everypath,tiltedformatedpaths{i},'IgnoreCase',true);
end
tiltedIndex=sum(Index,2);

% check non matches 
% C = setdiff(controlformatedpaths,everypath(logical(controlIndex)))

% C = setdiff(tiltedformatedpaths,everypath(logical(tiltedIndex)))
[C,I] = setdiff(everypath(logical(tiltedIndex)),tiltedformatedpaths);

for i=1:length(C)
   negindex(:,i)=contains(everypath,C{i},'IgnoreCase',true);
end
negindex=sum(negindex,2);
tiltedIndex=tiltedIndex-negindex;

clearvars -except controlIndex tiltedIndex ratIDs everypath

%% EXTRACT DATA [peak rate, coherence, information content, average rate]
data=importdata('/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/RH Place Cell Data Worksheet.4.xlsx');

All_data=[0 0 0 0 0];
for i=1:length(ratIDs)
    tempdat=data.data.(ratIDs{i});
    % extract waveform durations and fill in waveform duration for every trial
    waveform=tempdat(:,end);
    waveform(2:5:end)=waveform(1:5:end);
    waveform(3:5:end)=waveform(1:5:end);
    waveform(4:5:end)=waveform(1:5:end);
    waveform(5:5:end)=waveform(1:5:end);

    %format waveform duration into microseconds
    waveform=(waveform/33)*31.25;
    
    % store all data
    All_data=[All_data;[tempdat(:,5:8),waveform]];
end
All_data(1,:)=[];

% none of these measure should be negative
All_data=abs(All_data);

% coherence should be 0 to 1
All_data(find(All_data(:,2)>10),2)=All_data(All_data(:,2)>10,2)/100;
All_data(find(All_data(:,2)>1),2)=All_data(All_data(:,2)>1,2)/10;

%% PULL OUT CONTROL & TILTED 
% individual ids
control={'SR019','SR020','RH040'};
tilted={'SR011','SR012','RH048','RH026','SR021','SR024','SR023','SK107'};

% look through paths for ids
for i=1:length(control)
    index(:,i)=contains(everypath,control{i});
end
allcontrolidx=sum(index,2);

clear index

for i=1:length(tilted)
    index(:,i)=contains(everypath,tilted{i});
end
alltiltedidx=sum(index,2);

% pack into 3d matrix so that trial number is the 3 dim
allcontrol3d(:,:,1)=allcontrolidx(1:5:end,:);
allcontrol3d(:,:,2)=allcontrolidx(2:5:end,:);
allcontrol3d(:,:,3)=allcontrolidx(3:5:end,:);
allcontrol3d(:,:,4)=allcontrolidx(4:5:end,:);
allcontrol3d(:,:,5)=allcontrolidx(5:5:end,:);

alltilted3d(:,:,1)=alltiltedidx(1:5:end,:);
alltilted3d(:,:,2)=alltiltedidx(2:5:end,:);
alltilted3d(:,:,3)=alltiltedidx(3:5:end,:);
alltilted3d(:,:,4)=alltiltedidx(4:5:end,:);
alltilted3d(:,:,5)=alltiltedidx(5:5:end,:);

% same for actual data values
All_data3d(:,:,1)=All_data(1:5:end,:);
All_data3d(:,:,2)=All_data(2:5:end,:);
All_data3d(:,:,3)=All_data(3:5:end,:);
All_data3d(:,:,4)=All_data(4:5:end,:);
All_data3d(:,:,5)=All_data(5:5:end,:);

% pack data in 3d based on group
controldata=All_data3d(allcontrol3d(:,:,1)==1,:,:);
tilteddata=All_data3d(alltilted3d(:,:,1)==1,:,:);

% pack data in 3d based on group and if it is a place cell

controlIndex3d(:,:,1)=controlIndex(1:5:end,:);
controlIndex3d(:,:,2)=controlIndex(2:5:end,:);
controlIndex3d(:,:,3)=controlIndex(3:5:end,:);
controlIndex3d(:,:,4)=controlIndex(4:5:end,:);
controlIndex3d(:,:,5)=controlIndex(5:5:end,:);

tiltedIndex3d(:,:,1)=tiltedIndex(1:5:end,:);
tiltedIndex3d(:,:,2)=tiltedIndex(2:5:end,:);
tiltedIndex3d(:,:,3)=tiltedIndex(3:5:end,:);
tiltedIndex3d(:,:,4)=tiltedIndex(4:5:end,:);
tiltedIndex3d(:,:,5)=tiltedIndex(5:5:end,:);

controlplacecelldata=All_data3d(controlIndex3d(:,:,1)==1,:,:);
tiltedplacecelldata=All_data3d(tiltedIndex3d(:,:,1)==1,:,:);


everypath3d(:,:,1)=everypath(1:5:end,:);
everypath3d(:,:,2)=everypath(2:5:end,:);
everypath3d(:,:,3)=everypath(3:5:end,:);
everypath3d(:,:,4)=everypath(4:5:end,:);
everypath3d(:,:,5)=everypath(5:5:end,:);

allcontrolpaths=everypath3d(allcontrol3d(:,:,1)==1,:,:);
alltiltedpaths=everypath3d(alltilted3d(:,:,1)==1,:,:);

controlplacecellpaths=everypath3d(controlIndex3d(:,:,1)==1,:,1);
tiltedplacecellpaths=everypath3d(tiltedIndex3d(:,:,1)==1,:,1);



clearvars -except controlplacecelldata tiltedplacecelldata controldata tilteddata All_data3d everypath3d allcontrolpaths alltiltedpaths controlplacecellpaths tiltedplacecellpaths

save('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/RHPlaceCellDataWorksheet_Compiled.mat')

function out=formatpaths(id)
% convert names to same format created above to get a proper cross reference
for i=1:length(id)
    splitpath=strsplit(id{i},'/');
    splitdate=strsplit(splitpath{6},'_');
    splitname=strsplit(splitdate{2},' ');
    splitdateevenmore=strsplit(splitdate{1},'-');
    splitdateevenmore{3} = regexprep(splitdateevenmore{3},'^0*','');
    splitTT=strsplit(splitpath{8},'-');
    tempTT=strsplit(splitTT{1},' ');
    splitTT{1}=tempTT{1};
    cellnum=strsplit(splitTT{2},'.');
    cellnum=strsplit(cellnum{1},'_');
    cellnum{end} = regexprep(cellnum{end},'^0*','');
    formatedpaths{i,1}=[splitdateevenmore{1},' ',splitdateevenmore{2},'-',splitdateevenmore{3},' ',splitname{2},'/','ASCII','/',splitTT{1},' ','SS_',cellnum{end}];
end
out=upper(formatedpaths);
end