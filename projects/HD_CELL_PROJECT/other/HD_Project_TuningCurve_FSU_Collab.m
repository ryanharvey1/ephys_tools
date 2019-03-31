% HD_Project_TuningCurve_FSU_Collab
% opens .r text files containing HD cell data and creates tuning curves
%
%
% Ryan E Harvey 2018
%
% cd to data and get file names
clear;clc;close all
addpath(genpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/CircStat2012a'))
addpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis')

% path='/Users/ryanharvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis/HD_CELL_PROJECT/PHC HD Cells';
path='/Users/ryanharvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis/HD_CELL_PROJECT/other/PHC HD Cells';
cd(path)
filenames=dir;
filenames={filenames.name}';
filenames(contains(filenames,'._'))=[];
filenames(~contains(filenames,'.txt'))=[];


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
    
    % check firing in S2-S4
    if sum(sum(frames(:,7:9)))>0
        disp([filenames{i},' Has spikes in S2-S4'])
    end
    
%     frames(:,1)=linspace(0,(length(frames)/60),length(frames));
    

     frames(:,2:5)=frames(:,2:5)+1;
    
    % remove non-detects in the first 4 xy columns
    frames(frames(:,2)==255 | frames(:,3)==255 | frames(:,4)==255 | frames(:,5)==255,:)=NaN;
    
    [frames(:,2),frames(:,3)]=FixPos(frames(:,2),frames(:,3),frames(:,1),round(0.1667*60));

%     [thetaout]=smoothangle(frames(:,10),10);

    % expand frames to create spike binary
    spksover60hz=frames(frames(:,6)>1,:);
    
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
    
    clear tempframes EXP spksover60hz ii I
    
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
    
    %% COMPUTE STATS
    rlength = circ_r(angBins',hdTuning',deg2rad(6));
    
    %% save data
    id=regexprep(filenames{i},'[.txt]','');
    data.(id).frames_w_spk=framesEXP;
    data.(id).frames=frames;
    data.(id).tuningcurve=hdTuning;
    data.(id).rlength=rlength;
    

    %% PLOT
    tuningfig=figure; tuningfig.Color=[1 1 1];
    p=plot(hdTuning,'k')
    xlabel('Head Angle')
    ylabel('Firing Rate (hz)')
    title(['rlength: ',num2str(rlength)])
    set(p,'LineWidth',3)
    set(gca,'box','off','LineWidth',2,'XTick',linspace(0,60,7),'XTickLabel',linspace(0,60,7)*6,'FontSize',20,'FontWeight','bold')
    
    %print(tuningfig,'-dpng', '-r150',['/Users/ryanharvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis/HD_CELL_PROJECT/PHC HD Cells/TuningCurves/',id,'.png'])

 
    %%
end










