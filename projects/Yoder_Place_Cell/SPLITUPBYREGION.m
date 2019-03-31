% SPLITUPBYREGION
% Comment 2: p. 20, line 710 "We ?made no attempt to separate the
% recordings of neurons in area CA1 from those in dentate gyrus."
% While it seems unlikely that many cells were recorded in the dentate gyrus,
% it is difficult to compare across groups if the reader is completely unaware
% of the within group proportions of cells in CA1 vs DG. Is it possible any of
% these effects are actually due to differences in sampling in CA1vs DG across
% tilted vs control mice? Similarly, p. 18, line 570 seems to suggest that
% recordings might also include CA3, was CA3 also included and if so what were
% the relative proportions of cells from each hippocampal region in each group?

clear;clc;close all
addpath ('/Users/ryanharvey/GoogleDrive/MatlabDir/CircStat2012a','/Users/ryanharvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis')
load('/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/data5_SPLITUPBYREGION');
FigureLocation='/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/Figures_For_Paper';

reidx.sessXregionC=sessXregionC;
reidx.sessXregionT=sessXregionT;


%___________
% method 1: Ryan pipeline results 
% method 2: Jenny pipeline results
method=1;
%-----------
if method==1
    DispVarNames={'PeakRate','AVGRate','PercentActiveBins20','PercentActiveBins50','FieldWidth20','FieldWidth50','Sparsity','InfoContent','Coherence','FieldtoWall20per','FieldtoWall50per','PeakBin2Wall','BorderScore','fieldwidth','infieldFR','outfieldFR','Eccentricity','Velocity'};
    [Results]=analyze(sessionC,sessXregionC,resultsC,reidx,method);
    ca1c=Results.control.ca1;
    ca3c=Results.control.ca3;
    ca1t=Results.tilted.ca1;
    ca3t=Results.tilted.ca3;
    
    idx=zeros(1,length(DispVarNames));idx(1,[4,5,6,11])=1;
    DispVarNames(logical(idx))=[];
    ca1c(:,logical(idx),:)=[];
    ca3c(:,logical(idx),:)=[];
    ca1t(:,logical(idx),:)=[];
    ca3t(:,logical(idx),:)=[];
    
else
    DispVarNames={'PeakRate','nSpikes','OverallFR','NumbActiveBins','sparsity','InformationContent','Coherence','Field2Wall','borderScore','FieldWidth','infieldFR','outfieldFR','E','c','p'};
    [ResultsC]=analyze(sessionC,sessXregionC,resultsC,reidx,method);
    [ResultsT]=analyze(sessionT,sessXregionT,resultsT,reidx,method);
    ca1c=ResultsC.ca1;
    ca3c=ResultsC.ca3;
    ca1t=ResultsT.ca1;
    ca3t=ResultsT.ca3;
end

[ CvC ] = ScatterBox(ca1c(:,:,1),ca3c(:,:,1),{'Con CA1','Con CA3'},DispVarNames,1)
[ TvT ] = ScatterBox(ca1t(:,:,1),ca3t(:,:,1),{'Tilt CA1','Tilt CA3'},DispVarNames,1)

[ CvTca1 ] = ScatterBox(ca1c(:,:,1),ca1t(:,:,1),{'Con CA1','Tilt CA1'},DispVarNames,1)
[ CvTca3 ] = ScatterBox(ca3c(:,:,1),ca3t(:,:,1),{'Con CA3','Tilt CA3'},DispVarNames,1)

[ CvT ] = ScatterBox([ca1c(:,:,1);ca3c(:,:,1)],[ca1t(:,:,1);ca3t(:,:,1)],{'Con','Tilt'},DispVarNames,1)


% [ CvTca1vca3 ] = ScatterBox(ResultsC.ca1(:,:,1),ResultsT.=ca3(:,:,1),{'Con CA1','Tilt CA3'},DispVarNames,1)

ca1C=size(ca1c(:,:,1),1);
ca3C=size(ca3c(:,:,1),1);
ca1T=size(ca1t(:,:,1),1);
ca3T=size(ca3t(:,:,1),1);
[h,p, chi2stat,df] = prop_test([ca1C ca1T] , [ca1C+ca3C ca1T+ca3T], 0)
[h,p, chi2stat,df] = prop_test([ca3C ca3T] , [ca1C+ca3C ca1T+ca3T], 0)

disp([num2str(ca1C/(ca1C+ca3C)),' CA1 Control'])
disp([num2str(ca1T/(ca1T+ca3T)),' CA1 Tilted'])

% p = sampsizepwr('t2',[1.44 std(ca1c(:,8,1))],2.09,[],ca1C+ca3C)
% p = sampsizepwr('t2',[nanmean(ca1t(:,8,1)) nanstd(ca1t(:,8,1))],nanmean(ca3t(:,8,1)),[],ca1T+ca3T)
for i=7:9
    disp('control')
    disp([num2str(round(mean(ca3c(:,i,1)),2)),' ± ',num2str(round(std(ca3c(:,i,1))/sqrt(length(ca3c(:,i,1))),2)),'    ',...
        num2str(round(mean(ca3c(:,i,2)),2)),' ± ',num2str(round(std(ca3c(:,i,2))/sqrt(length(ca3c(:,i,2))),2)),'    ',...
        num2str(round(mean(ca3c(:,i,3)),2)),' ± ',num2str(round(std(ca3c(:,i,3))/sqrt(length(ca3c(:,i,3))),2)),'    ',...
        num2str(round(mean(ca3c(:,i,4)),2)),' ± ',num2str(round(std(ca3c(:,i,4))/sqrt(length(ca3c(:,i,4))),2)),'    ',...
        num2str(round(mean(ca3c(:,i,5)),2)),' ± ',num2str(round(std(ca3c(:,i,5))/sqrt(length(ca3c(:,i,5))),2))])
    disp('tilted')
    disp([num2str(round(mean(ca3t(:,i,1)),2)),' ± ',num2str(round(std(ca3t(:,i,1))/sqrt(length(ca3t(:,i,1))),2)),'    ',...
        num2str(round(mean(ca3t(:,i,2)),2)),' ± ',num2str(round(std(ca3t(:,i,2))/sqrt(length(ca3t(:,i,2))),2)),'    ',...
        num2str(round(mean(ca3t(:,i,3)),2)),' ± ',num2str(round(std(ca3t(:,i,3))/sqrt(length(ca3t(:,i,3))),2)),'    ',...
        num2str(round(mean(ca3t(:,i,4)),2)),' ± ',num2str(round(std(ca3t(:,i,4))/sqrt(length(ca3t(:,i,4))),2)),'    ',...
        num2str(round(mean(ca3t(:,i,5)),2)),' ± ',num2str(round(std(ca3t(:,i,5))/sqrt(length(ca3t(:,i,5))),2))])
end
% napprox = sampsizepwr('p',nanmean(ca1t(:,8,1)),nanmean(ca3t(:,8,1)),0.8)

function [results]=analyze(sess,sessXregion,data,reidx,method)
if method==2
    sessXregion=sessXregion(~cellfun('isempty',sessXregion(:,1)),:);
    sessXregion=sessXregion(~cellfun('isempty',sessXregion(:,2)),:);
    
    sess=sess(~cellfun('isempty',sess(:,1)),:);
    
    sessXregion(ismember(sessXregion(:,2),'dentate'),2)={'CA3'};
    
    CA1=sessXregion(ismember(sessXregion(:,2),'CA1'),1);
    CA3=sessXregion(ismember(sessXregion(:,2),'CA3'),1);
    
    idxca1=ismember(sess(:),CA1);
    idxca3=ismember(sess(:),CA3);
    
    results.ca1=data(idxca1,:,:);
    results.ca3=data(idxca3,:,:);
elseif method==1
    load('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/GW_Field_StatsPlaceCells_Tilted_Mice.mat')
    load('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/NewData.mat')
    load('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/NewData2.mat');
    newdata2=ans;
    allcells=newdata.textdata(2:end,1);
    placecells=data.textdata.Field_StatsPlaceCells0x2DTilted(2:end,1);
    for i=1:length(placecells);newmeasures(i,:)=newdata.data(strcmpi(placecells(i),allcells),13:16);end
    allcells2=newdata2.textdata(2:end,1);for i=1:length(placecells);newmeasures2(i,:)=newdata2.data(strcmpi(placecells(i),allcells2),17:18);end
    data.data.Field_StatsPlaceCells0x2DTilted=[data.data.Field_StatsPlaceCells0x2DTilted,newmeasures,newmeasures2];
    control=data.data.Field_StatsPlaceCells0x2DTilted(ismember(data.textdata.Field_StatsPlaceCells0x2DTilted(2:end,2),'Control'),:);
    tilted=data.data.Field_StatsPlaceCells0x2DTilted(ismember(data.textdata.Field_StatsPlaceCells0x2DTilted(2:end,2),'Tilted'),:);
    for i=1:5;Group1(:,:,i)=control(i:5:end,:);Group2(:,:,i)=tilted(i:5:end,:);end
    % ADD ECCENTRCITY
    load('/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/New_eccentricity.mat')
    [r,c,d]=size(Group1);for i=1:5;Group1(:,c+1,i)=eccentricity.ec(:,1,i);end
    [r,c,d]=size(Group2);for i=1:5;Group2(:,c+1,i)=eccentricity.et(:,1,i);end
    % ADD Velocity
    load('/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/VelocityOverSessions.mat')
    [r,c,d]=size(Group1);for i=1:5;Group1(:,c+1,i)=VelocityOverSessions.controlresult(:,1,i);end
    [r,c,d]=size(Group2);for i=1:5;Group2(:,c+1,i)=VelocityOverSessions.tiltedresult(:,1,i);end
    VarNames=[data.textdata.Field_StatsPlaceCells0x2DTilted(1,3:end),newdata.textdata(1,end-3:end),newdata2.textdata(1,end-1:end)];
    VarNames=regexprep(VarNames,'_','','emptymatch');VarNames=[VarNames,'Eccentricity','Velocity'];controldisp=Group1(:,15:16,2); Group1(:,15:16,:)=[];tiltdisp=Group2(:,15:16,2); Group2(:,15:16,:)=[];DispVarNames=VarNames(15:16); VarNames(15:16)=[];
    
    
    pathC=data.textdata.Field_StatsPlaceCells0x2DTilted(ismember(data.textdata.Field_StatsPlaceCells0x2DTilted(:,2),'Control'),1);
    pathT=data.textdata.Field_StatsPlaceCells0x2DTilted(ismember(data.textdata.Field_StatsPlaceCells0x2DTilted(:,2),'Tilted'),1);
    pathC=pathC([1:5:length(pathC)],1);
    pathT=pathT([1:5:length(pathT)],1);
    
    for i=1:length(pathC)
        temppath=strsplit(pathC{i},filesep);
        pathC(i)=temppath(6);
    end
    for i=1:length(pathT)
        temppath=strsplit(pathT{i},filesep);
        pathT(i)=temppath(6);
    end
    
    sessXregion=reidx.sessXregionC;
    
    sessXregion=sessXregion(~cellfun('isempty',sessXregion(:,1)),:);
    sessXregion=sessXregion(~cellfun('isempty',sessXregion(:,2)),:);
    sessXregion(ismember(sessXregion(:,2),'dentate'),2)={'CA3'};
    
    CA1=sessXregion(ismember(sessXregion(:,2),'CA1'),1);
    CA3=sessXregion(ismember(sessXregion(:,2),'CA3'),1);
    
    idxca1=ismember(pathC(:),CA1);
    idxca3=ismember(pathC(:),CA3);
    
    results.control.ca1=Group1(idxca1,:,:);
    results.control.ca3=Group1(idxca3,:,:);
    
    %/////////////////////////////////////////////////////////////////////
    
    sessXregion=reidx.sessXregionT;
    
    sessXregion=sessXregion(~cellfun('isempty',sessXregion(:,1)),:);
    sessXregion=sessXregion(~cellfun('isempty',sessXregion(:,2)),:);
    sessXregion(ismember(sessXregion(:,2),'dentate'),2)={'CA3'};
    
    CA1=sessXregion(ismember(sessXregion(:,2),'CA1'),1);
    CA3=sessXregion(ismember(sessXregion(:,2),'CA3'),1);
    
    idxca1=ismember(pathT(:),CA1);
    idxca3=ismember(pathT(:),CA3);
    
    results.tilted.ca1=Group2(idxca1,:,:);
    results.tilted.ca3=Group2(idxca3,:,:);
end
end


