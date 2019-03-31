% HD_Project_TuningCurve_v2
% opens .r text files containing HD cell data and creates tuning curves
%
%
% Ryan E Harvey 2018
%
% cd to data and get file names
clear;clc;close all
com=which('HD_Project_TuningCurve_v2');
com=strsplit(com,filesep);
basedir=[com{1},filesep,'Users',filesep,com{3},filesep,'GoogleDrive',filesep,'MatlabDir'];
addpath(genpath([basedir,filesep,filesep,'CircStat2012a']),...
    genpath([basedir,filesep,'BClarkToolbox',filesep, 'Analysis']));

addpath(genpath('C:\Users\Ben Clark''s Lab\Google Drive\MatlabDir\BClarkToolbox\Analysis'));
addpath(genpath('C:\Users\Ben Clark''s Lab\Google Drive\MatlabDir\CircStat2012a'));


path='d:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\Ben_HDCProject\Data';


cd(path)
files=dir( '**/*.r');
folders=unique({files.folder});

folders(contains(folders,'Timestamps'))=[];
folders(contains(folders,'timestamp'))=[];


for i=1:length(folders)
    area=strsplit(folders{i},filesep);
    area=area{end};
    area=erase(area,'_SpikePositionFiles');
    areas{i}=erase(area,'_');
end
clear area files path

for a=1:length(folders)
    cd(folders{a})
    filenames=dir('*.r');
    filenames={filenames.name}';
    filenames(contains(filenames,'._'))=[];
    
    for i=1:length(filenames)
        % OPEN TEXT FILE AND CLEAN DATA
        disp(['Running:  ',filenames{i}])
        fileID=fopen(filenames{i},'r');
        dataArray=textscan(fileID,'%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]','Delimiter',...
            '\t','TextType','string','EmptyValue',NaN,'HeaderLines',2-1,'ReturnOnError',false,'EndOfLine','\r\n');
        fclose(fileID);
        frames=[dataArray{1:end-1}];
        clear dataArray fileID
        % remove last column (all NaNs)
        frames(:,end)=[];
        
        % convert ts to secs
        samplerate=60;
        frames(:,1)=linspace(0,(length(frames(:,1))/samplerate),length(frames(:,1)));
        
        truesamplerateTS=frames(:,1);
        
        % remove non-detects in the first 4 xy columns
        frames(frames(:,2)==0 | frames(:,3)==0 | frames(:,4)==0 | frames(:,5)==0,:)=[];
        
        % expand frames to create spike binary
        spksover60hz=frames(frames(:,6)>1,:);
        
        idx=find(frames(:,6)>1);
        ts=[];
        for ii=1:size(spksover60hz,1)
            if idx(ii)+1>length(frames)
                tstemp=linspace(frames(idx(ii)-1,1),frames(idx(ii),1),frames(idx(ii),6)+2);
            elseif idx(ii)-1<length(frames)
                tstemp=linspace(frames(idx(ii),1),frames(idx(ii)+1,1),frames(idx(ii),6)+2);
            else
                tstemp=linspace(frames(idx(ii)-1,1),frames(idx(ii)+1,1),frames(idx(ii),6)+2);
            end
            tstemp(1)=[];
            tstemp(end)=[];
            ts=[ts;tstemp'];
        end
        clear tstemp idx
        ts=sort([frames(frames(:,6)<=1,1);ts]);
        
        % remove frames with > spike binary 1 for later recombining
        tempframes=frames;
        tempframes(tempframes(:,6)>1,:)=[];
        
        EXP=zeros(1,size(spksover60hz,2));
        for ii=1:size(spksover60hz,1)
            EXP=[EXP;repmat(spksover60hz(ii,:),spksover60hz(ii,6),1)];
        end
        EXP(1,:)=[];
        EXP(:,6)=ones(size(EXP,1),1);
        
        framesEXP=[tempframes;EXP];
        
        [~,I]=sort(framesEXP(:,1));
        
        framesEXP=framesEXP(I,:);
        
        framesEXP(:,1)=ts;
        
        HDdata.(areas{a}).frames{i,1}=framesEXP;
        
        clear tempframes EXP spksover60hz ii I
        
        %% ////////////// HD CELL ANAYSIS ///////////////////
        %
        % STABILITY
        [HDdata.(areas{a}).within_Coeff(i,1),~,within]=within_HDstability(framesEXP,framesEXP(framesEXP(:,6)==0,:),60,10,6);
        
        HDdata.(areas{a}).withinTuningCurves1(i,:)=within(1,:);
        HDdata.(areas{a}).withinTuningCurves2(i,:)=within(2,:);
        HDdata.(areas{a}).withinTuningCurves3(i,:)=within(3,:);
        HDdata.(areas{a}).withinTuningCurves4(i,:)=within(4,:);
        
        
        %% CREATE TUNING CURVE
        % 6 degree bins
        da=pi/30;
        angBins=[da/2:da:2*pi-da/2];
        % Occupancy
        histAng=hist(frames(:,10),angBins);
        % Number of spikes per bin
        spkPerAng=hist(framesEXP(framesEXP(:,6)==1,10),angBins);
        % Tuning
        hdTuning=(spkPerAng./histAng)*60;
        % remove nan & inf
        hdTuning(isnan(hdTuning) | isinf(hdTuning))=0;
        
        clear da histAng spkPerAng
        
        %% COMPUTE R LENGTH
        rlength(i,1)=circ_r(angBins',hdTuning',deg2rad(6));
        
        %% CIRCSHIFT SO PEAK IS IN THE MIDDLE
        %         [M,I]=max(hdTuning);
        %         middlebin=round(median(1:length(hdTuning)));
        %
        %         hdTuningshift=circshift(hdTuning,(middlebin-I)-1);
        
        HDdata.(areas{a}).rawTuningCurve(i,:)=hdTuning;
        
        clear M I middlebin
        
        %% Normalize & Collect tuning curve
        ALL_hdTuning(i,:)=rescale(hdTuning,0,1);
        %                 ALL_hdTuning(i,:)=rescale(hdTuningshift,0,1);
        
        % Save Temp File
        %       save(['D:\Projects\HDtempData',filesep,areas{a},'_iteration_',num2str(i),'.mat'],'HDdata','-v7.3')
    end
    %% SORT AND SAVE
    [~,pdIdx]=max(ALL_hdTuning,[],2);
    [~,I]=sort(pdIdx);
    ALL_hdTuningsorted=ALL_hdTuning(I,:);
    
    HDdata.(areas{a}).TuningCurve=ALL_hdTuningsorted;
    HDdata.(areas{a}).TuningIdx=I;
    HDdata.(areas{a}).RLength=rlength;
    
    clear ALL_hdTuning rlength
end


%% 4-Quarter Stability Figure

ATN=[HDdata.ATN.rawTuningCurve];
PoS=[HDdata.PoSDeep.rawTuningCurve;HDdata.PoSSup.rawTuningCurve];
MEC=[HDdata.MECDeepLayers.rawTuningCurve;HDdata.MECSupLayers.rawTuningCurve];
PaS=[HDdata.PaSDeep.rawTuningCurve;HDdata.PaSSup.rawTuningCurve];

ATN1=[HDdata.ATN.withinTuningCurves1];
PoS1=[HDdata.PoSDeep.withinTuningCurves1;HDdata.PoSSup.withinTuningCurves1];
MEC1=[HDdata.MECDeepLayers.withinTuningCurves1;HDdata.MECSupLayers.withinTuningCurves1];
PaS1=[HDdata.PaSDeep.withinTuningCurves1;HDdata.PaSSup.withinTuningCurves1];

ATN2=[HDdata.ATN.withinTuningCurves2];
PoS2=[HDdata.PoSDeep.withinTuningCurves2;HDdata.PoSSup.withinTuningCurves2];
MEC2=[HDdata.MECDeepLayers.withinTuningCurves2;HDdata.MECSupLayers.withinTuningCurves2];
PaS2=[HDdata.PaSDeep.withinTuningCurves2;HDdata.PaSSup.withinTuningCurves2];

ATN3=[HDdata.ATN.withinTuningCurves3];
PoS3=[HDdata.PoSDeep.withinTuningCurves3;HDdata.PoSSup.withinTuningCurves3];
MEC3=[HDdata.MECDeepLayers.withinTuningCurves3;HDdata.MECSupLayers.withinTuningCurves3];
PaS3=[HDdata.PaSDeep.withinTuningCurves3;HDdata.PaSSup.withinTuningCurves3];

ATN4=[HDdata.ATN.withinTuningCurves4];
PoS4=[HDdata.PoSDeep.withinTuningCurves4;HDdata.PoSSup.withinTuningCurves4];
MEC4=[HDdata.MECDeepLayers.withinTuningCurves4;HDdata.MECSupLayers.withinTuningCurves4];
PaS4=[HDdata.PaSDeep.withinTuningCurves4;HDdata.PaSSup.withinTuningCurves4];

%Get first quarter index
[~,pdIdx]=max(ATN,[],2);
[~,ATNidx]=sort(pdIdx);
ATN_sorted=ATN(ATNidx,:);
for i=1:size(ATN_sorted,1)
    plot_Overall_ATN(i,:)=rescale(ATN_sorted(i,:),0,1);
end

[~,pdIdx]=max(PoS,[],2);
[~,PoSidx]=sort(pdIdx);
PoS_sorted=PoS(PoSidx,:);
for i=1:size(PoS_sorted,1)
    plot_Overall_PoS(i,:)=rescale(PoS_sorted(i,:),0,1);
end

[~,pdIdx]=max(MEC,[],2);
[~,MECidx]=sort(pdIdx);
MEC_sorted=MEC(MECidx,:);
for i=1:size(MEC_sorted,1)
    plot_Overall_MEC(i,:)=rescale(MEC_sorted(i,:),0,1);
end

[~,pdIdx]=max(PaS,[],2);
[~,PaSidx]=sort(pdIdx);
PaS_sorted=PaS(PaSidx,:);
for i=1:size(PaS_sorted,1)
    plot_Overall_PaS(i,:)=rescale(PaS_sorted(i,:),0,1);
end

[~,pdIdx]=max(ATN1,[],2);
[~,ATN1idx]=sort(pdIdx);

[~,pdIdx]=max(PoS1,[],2);
[~,PoS1idx]=sort(pdIdx);

[~,pdIdx]=max(MEC1,[],2);
[~,MEC1idx]=sort(pdIdx);

[~,pdIdx]=max(PaS1,[],2);
[~,PaS1idx]=sort(pdIdx);

%create matrices
ATNQtr1Tuning_sorted=ATN1(ATN1idx,:); ATNQ1_corr=corr(ATNQtr1Tuning_sorted,ATNQtr1Tuning_sorted);
PoSQtr1Tuning_sorted=PoS1(PoS1idx,:); PoSQ1_corr=corr(PoSQtr1Tuning_sorted,PoSQtr1Tuning_sorted);
MECQtr1Tuning_sorted=MEC1(MEC1idx,:); MECQ1_corr=corr(MECQtr1Tuning_sorted,MECQtr1Tuning_sorted);
PaSQtr1Tuning_sorted=PaS1(PaS1idx,:); PaSQ1_corr=corr(PaSQtr1Tuning_sorted,PaSQtr1Tuning_sorted);


ATNQtr2Tuning_sorted=ATN2(ATN1idx,:); ATNQ2_corr=corr(ATNQtr2Tuning_sorted,ATNQtr2Tuning_sorted);
PoSQtr2Tuning_sorted=PoS2(PoS1idx,:); PoSQ2_corr=corr(PoSQtr2Tuning_sorted,PoSQtr2Tuning_sorted);
MECQtr2Tuning_sorted=MEC2(MEC1idx,:); MECQ2_corr=corr(MECQtr2Tuning_sorted,MECQtr2Tuning_sorted);
PaSQtr2Tuning_sorted=PaS2(PaS1idx,:); PaSQ2_corr=corr(PaSQtr2Tuning_sorted,PaSQtr2Tuning_sorted);


ATNQtr3Tuning_sorted=ATN3(ATN1idx,:); ATNQ3_corr=corr(ATNQtr3Tuning_sorted,ATNQtr3Tuning_sorted);
PoSQtr3Tuning_sorted=PoS3(PoS1idx,:); PoSQ3_corr=corr(PoSQtr3Tuning_sorted,PoSQtr3Tuning_sorted);
MECQtr3Tuning_sorted=MEC3(MEC1idx,:); MECQ3_corr=corr(MECQtr3Tuning_sorted,MECQtr3Tuning_sorted);
PaSQtr3Tuning_sorted=PaS3(PaS1idx,:); PaSQ3_corr=corr(PaSQtr3Tuning_sorted,PaSQtr3Tuning_sorted);


ATNQtr4Tuning_sorted=ATN4(ATN1idx,:); ATNQ4_corr=corr(ATNQtr4Tuning_sorted,ATNQtr4Tuning_sorted);
PoSQtr4Tuning_sorted=PoS4(PoS1idx,:); PoSQ4_corr=corr(PoSQtr4Tuning_sorted,PoSQtr4Tuning_sorted);
MECQtr4Tuning_sorted=MEC4(MEC1idx,:); MECQ4_corr=corr(MECQtr4Tuning_sorted,MECQtr4Tuning_sorted);
PaSQtr4Tuning_sorted=PaS4(PaS1idx,:); PaSQ4_corr=corr(PaSQtr4Tuning_sorted,PaSQtr4Tuning_sorted);


all_ATN_corr=corr(plot_Overall_ATN,plot_Overall_ATN);
all_PoS_corr=corr(plot_Overall_PoS,plot_Overall_PoS);
all_MEC_corr=corr(plot_Overall_MEC,plot_Overall_MEC);
all_PaS_corr=corr(plot_Overall_PaS,plot_Overall_PaS);


%combine for plot

plot_ATN=[plot_Overall_ATN nan(length(ATNQtr1Tuning_sorted),5)...
    ATNQtr1Tuning_sorted nan(length(ATNQtr1Tuning_sorted),1)...
    ATNQtr2Tuning_sorted nan(length(ATNQtr1Tuning_sorted),1)...
    ATNQtr3Tuning_sorted nan(length(ATNQtr1Tuning_sorted),1) ...
    ATNQtr4Tuning_sorted];
plot_PoS=[plot_Overall_PoS nan(length(PoSQtr1Tuning_sorted),5)...
    PoSQtr1Tuning_sorted nan(length(PoSQtr1Tuning_sorted),1)...
    PoSQtr2Tuning_sorted nan(length(PoSQtr1Tuning_sorted),1)...
    PoSQtr3Tuning_sorted nan(length(PoSQtr1Tuning_sorted),1)...
    PoSQtr4Tuning_sorted];
plot_PaS=[plot_Overall_PaS nan(length(PaSQtr1Tuning_sorted),5)...
    PaSQtr1Tuning_sorted nan(length(PaSQtr1Tuning_sorted),1)...
    PaSQtr2Tuning_sorted nan(length(PaSQtr1Tuning_sorted),1)...
    PaSQtr3Tuning_sorted nan(length(PaSQtr1Tuning_sorted),1)...
    PaSQtr4Tuning_sorted];
plot_MEC=[plot_Overall_MEC nan(length(MECQtr1Tuning_sorted),5)...
    MECQtr1Tuning_sorted nan(length(MECQtr1Tuning_sorted),1)...
    MECQtr2Tuning_sorted nan(length(MECQtr1Tuning_sorted),1)...
    MECQtr3Tuning_sorted nan(length(MECQtr1Tuning_sorted),1)...
    MECQtr4Tuning_sorted];

%POP VEC
fig=figure; fig.Color=[1 1 1];
subplot(4,1,1)
imAlpha=ones(size(plot_ATN));
imAlpha(isnan(plot_ATN))=0;
imagesc(plot_ATN,'AlphaData',imAlpha); colormap parula(255);box off;
ylabel('Cells')
set(gca,'XTick',[30.7,95.7,155.7,225.7, 280.7],'XTickLabel',[],...
    'YTickLabel',[],'FontSize',20,'FontWeight','bold','TickLength',[0;0])
subplot(4,1,2)
imAlpha=ones(size(plot_PoS));
imAlpha(isnan(plot_PoS))=0;
imagesc(plot_PoS,'AlphaData',imAlpha); colormap parula(255);box off;
ylabel('Cells')
set(gca,'XTick',[30.7,95.7,155.7,225.7,280.7],'XTickLabel',[],...
    'YTickLabel',[],'FontSize',20,'FontWeight','bold','TickLength',[0;0])
subplot(4,1,3)
imAlpha=ones(size(plot_PaS));
imAlpha(isnan(plot_PaS))=0;
imagesc(plot_PaS,'AlphaData',imAlpha); colormap parula(255);box off;
ylabel('Cells')
set(gca,'XTick',[30.7,95.7,155.7,225.7, 280.7],'XTickLabel',[],...
    'YTickLabel',[],'FontSize',20,'FontWeight','bold','TickLength',[0;0])
subplot(4,1,4)
imAlpha=ones(size(plot_MEC));
imAlpha(isnan(plot_MEC))=0;
imagesc(plot_MEC,'AlphaData',imAlpha); colormap parula(255);box off;
ylabel('Cells')
set(gca,'XTick',[30.7,95.7,155.7,220.7, 280.7],'XTickLabel',{'Overall','First Quarter','Second Quarter','Third Quarter',...
    'Fourth Quarter'},'YTickLabel',[],'FontSize',20,'FontWeight','bold','TickLength',[0;0]);


%AUTOCORR PLOT
fig=figure; fig.Color=[1 1 1];
subplot(4,5,1)
imAlpha=ones(size(plot_ATN));
imAlpha(isnan(plot_ATN))=0;
imagesc(plot_ATN,'AlphaData',imAlpha); colormap parula(255);box off;
ylabel('Cells')
set(gca,'XTick',[30.7,95.7,155.7,225.7, 280.7],'XTickLabel',[],...
    'YTickLabel',[],'FontSize',20,'FontWeight','bold','TickLength',[0;0])
subplot(4,1,2)
imAlpha=ones(size(plot_PoS));
imAlpha(isnan(plot_PoS))=0;
imagesc(plot_PoS,'AlphaData',imAlpha); colormap parula(255);box off;
ylabel('Cells')
set(gca,'XTick',[30.7,95.7,155.7,225.7,280.7],'XTickLabel',[],...
    'YTickLabel',[],'FontSize',20,'FontWeight','bold','TickLength',[0;0])
subplot(4,1,3)
imAlpha=ones(size(plot_PaS));
imAlpha(isnan(plot_PaS))=0;
imagesc(plot_PaS,'AlphaData',imAlpha); colormap parula(255);box off;
ylabel('Cells')
set(gca,'XTick',[30.7,95.7,155.7,225.7, 280.7],'XTickLabel',[],...
    'YTickLabel',[],'FontSize',20,'FontWeight','bold','TickLength',[0;0])
subplot(4,1,4)
imAlpha=ones(size(plot_MEC));
imAlpha(isnan(plot_MEC))=0;
imagesc(plot_MEC,'AlphaData',imAlpha); colormap parula(255);box off;
ylabel('Cells')
set(gca,'XTick',[30.7,95.7,155.7,220.7, 280.7],'XTickLabel',{'Overall','First Quarter','Second Quarter','Third Quarter',...
    'Fourth Quarter'},'YTickLabel',[],'FontSize',20,'FontWeight','bold','TickLength',[0;0]);
subplot(4,2,1)
imagesc(plot_Overall_ATN); colormap parula(255);box off;
ylabel('Cells')
set(gca,'XTick',[30.7],'XTickLabel',[],...
    'YTickLabel',[],'FontSize',20,'FontWeight','bold','TickLength',[0;0])
subplot(4,2,3)
imagesc(plot_Overall_PoS); colormap parula(255);box off;
ylabel('Cells')
set(gca,'XTick',[30.7],'XTickLabel',[],...
    'YTickLabel',[],'FontSize',20,'FontWeight','bold','TickLength',[0;0])
subplot(4,2,5)
imagesc(plot_Overall_PaS); colormap parula(255);box off;
ylabel('Cells')
set(gca,'XTick',[30.7],'XTickLabel',[],...
    'YTickLabel',[],'FontSize',20,'FontWeight','bold','TickLength',[0;0])
subplot(4,2,7)
imagesc(plot_Overall_MEC); colormap parula(255);box off;
ylabel('Cells')
set(gca,'XTick',[30.7],'XTickLabel',{'Overall'},...
    'YTickLabel',[],'FontSize',20,'FontWeight','bold','TickLength',[0;0]);


