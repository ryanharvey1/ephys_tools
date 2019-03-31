% MWM_IRR
% clear;close all;clc
load('C:\Users\Ben Clark''s Lab\Google Drive\MATLAB\MWM\oldData\NavDay1T3.mat')
% keep('pathResults','xposRot','yposRot')
ID=fieldnames(pathResults);
SEX=repmat([1 1 1 1 1 1 1 0 0 0 0 0 0 0 1 0 1 1 1 1 1 0 0 0 0 0 0]',20,1); %for T3
GENOTYPE=repmat([1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0]',20,1); %for T3
TRIAL=sortrows(repmat([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20]',27,1),1); %27 or 28=number of subjects

j=1;
% % NumofStrat=0:11;
for ii=1:length(ID)
trial=(strsplit(ID{ii,:},'_'));
strategytemp=pathResults.(ID{ii}).strategy;
    x=0; y=0; rad=75; th = 0:pi/179.5:2*pi;
    xunit=rad*cos(th)+x;yunit=rad*sin(th)+y;
    
    for i=1:length(strategytemp)
        fig=figure; fig.Color=[1 1 1];set(fig,'Position',[1 1 581 606]);
        plot(xunit,yunit,'k');hold on;axis image;box off;axis off;
        rectangle('Position',[xposRot yposRot 16 16]);hold on
        plot(pathResults.(ID{ii}).segments(i).segments(:,1),pathResults.(ID{ii}).segments(i).segments(:,2));
        scatter(pathResults.(ID{ii}).segments(i).segments(1,1),pathResults.(ID{ii}).segments(i).segments(1,2), 80,'filled', 'r')
%         strategy=(pathResults.(ID{ii}).strategy(i)+1);
        label(j,1)=ii;
        label(j,2)=i;
        an=input('what''s the strategy?    ');
        while ~ismember(an,1:12)
            disp('oops, must be a number between 1-12')
            an=input('what''s the strategy?    ');
        end
        while isempty(an)
            an=input('what''s the strategy?    ');
            if ~ismember(an,1:12)
                    disp('oops, must be a number between 1-12')
                    an=input('what''s the strategy?    ');
            end
            
       end
        label(j,3)=an;
        if an==12; ml=input('what are the multiple strategies?    ');
            multi.(['multi', num2str(j)])=ml;
        end
        j=j+1;
        if ii >=length(ID)*.25 && ii <=length(ID)*.26;disp('You have labeled 25% of data');end
        if ii >=length(ID)*.50 && ii <=length(ID)*.51;disp('You have labeled 50% of data');end
        if ii >=length(ID)*.75 && ii <=length(ID)*.76;disp('You have labeled 75% of data');end
        if ii >=length(ID)*.95 && ii <=length(ID)*.96;disp('Almost there!');end

        close all
    end
end
% Compile Labeled Data
load('T3abels.mat'); load('T3multiProbe.mat'); load('ProbeID_T3')   %load applicable data files
% DEFINE GROUPS
TgM=[1:7,15]; TgF=[8:14,16]; WtM=17:22; WtF=23:28;

% GET SEGMENT ID
segID=[]; trialID=[]; Genotype=[]; Sex=[];
for ii=1:length(ID)
segID=[segID;repmat(ID(ii),size(pathResults.(ID{ii}).segmentStats,1),1)];
tempTrialID=str2double(regexp(ID{ii,:},'(?<=t\s*)\d*', 'match'));
trialID=[trialID; repmat(tempTrialID,size(pathResults.(ID{ii}).segmentStats,1),1)];
Genotype=[Genotype;repmat(GENOTYPE (ii),size(pathResults.(ID{ii}).segmentStats,1),1)];
Sex=[Sex; repmat(SEX(ii),size(pathResults.(ID{ii}).segmentStats,1),1)];
end

%COMPILE INDEX WITH LABELS
fullLabel=[trialID Genotype Sex label];   %concatenate group index with data

%COMPILE MULTI STRATEGIES 
name=fieldnames(multi);
for k=1:length(fieldnames(multi)) %cycle through fieldnames
    tempStrat=multi.(name{k});    % pull strategy data
    pullIndex=str2double(regexp(name{k,:},'(?<=i\s*)\d*', 'match'));  %pull index of strategy
    tempNewRow=fullLabel(pullIndex,:);
    for j=1:size(tempStrat,2)
        newRow(j,:)=[tempNewRow(:,1:end-1) tempStrat(:,j)];
    end
    fullLabel=[fullLabel; newRow];
end

fullLabel=sortrows(fullLabel,4);
fullLabel(fullLabel(:,6)==12,:)=[];


% COMPILE FREQ BY TRIAL AND SUBJECT FOR CONVOLUTION ANALYSIS COMPUTATIONS
 for l=1:length(unique(fullLabel(:,4)))
     temp=fullLabel(fullLabel(:,4)==l,6);
        for ii=1:11
            por(ii,1)=sum(temp(:,1)==ii)/size(temp,1);
        end
     allcollect(l,:)=por';
 end

% REMOVE ZEROS FROM ALLCOLLECT
 allcollect(allcollect==0)=NaN;
 
 % COMPILE PATH LENGTH
 pathLength=[];
 load('workspaceDay1T2','mwmTrainOutput'); pathLength=[pathLength; mwmTrainOutput{:,5}];
 load('workspaceDay2T2','mwmTrainOutput'); pathLength=[pathLength; mwmTrainOutput{:,5}];
 load('workspaceDay3T2','mwmTrainOutput'); pathLength=[pathLength; mwmTrainOutput{:,5}];
 load('workspaceDay4T2','mwmTrainOutput'); pathLength=[pathLength; mwmTrainOutput{:,5}];
 load('workspaceDay5T2','mwmTrainOutput'); pathLength=[pathLength; mwmTrainOutput{:,5}];
 pathLength=repmat(pathLength,1,11);
 
 % TgM  FREQ Total
    temp=fullLabel(fullLabel(:,2)==1 & fullLabel(:,3)==1,6);
    for ii=1:11
        tgmTotal(ii,1)=sum(temp(:,1)==ii)/size(temp,1);
    end
    clear temp;
% WtM  FREQ Total
    temp=fullLabel(fullLabel(:,2)==0 & fullLabel(:,3)==1,6);
    for ii=1:11
        wtmTotal(ii,1)=sum(temp(:,1)==ii)/size(temp,1);
    end
clear temp;
%  TgF FREQ Total
    temp=fullLabel(fullLabel(:,2)==1 & fullLabel(:,3)==0,6);
    for ii=1:11
        tgfTotal(ii,1)=sum(temp(:,1)==ii)/size(temp,1);
    end
clear temp;
%  WTF FREQ Total
    temp=fullLabel(fullLabel(:,2)==0 & fullLabel(:,3)==0,6);
    for ii=1:11
        wtfTotal(ii,1)=sum(temp(:,1)==ii)/size(temp,1);
    end
clear temp;
    T3freqTotal=[wtmTotal tgmTotal wtfTotal tgfTotal];
 
% TgM  FREQ BY TRIAL
for i=1:length(unique(trialID))
    temp=fullLabel(fullLabel(:,2)==1 & fullLabel(:,3)==1 & fullLabel(:,1)==i,6);
    for ii=1:11
        por(ii,1)=sum(temp(:,1)==ii)/size(temp,1);
    end
    tgmcollect(i,:)=por';
end
%
% WtM  FREQ BY TRIAL
for i=1:length(unique(trialID))
    temp=fullLabel(fullLabel(:,2)==0 & fullLabel(:,3)==1 & fullLabel(:,1)==i,6);
    for ii=1:11
        por(ii,1)=sum(temp(:,1)==ii)/size(temp,1);
    end
    wtmcollect(i,:)=por';
end

%  TgF FREQ BY TRIAL
for i=1:length(unique(trialID))
    temp=fullLabel(fullLabel(:,2)==1 & fullLabel(:,3)==0 & fullLabel(:,1)==i,6);
    for ii=1:11
        por(ii,1)=sum(temp(:,1)==ii)/size(temp,1);
    end
    tgfcollect(i,:)=por';
end

%  WTF FREQ BY TRIAL
for i=1:length(unique(trialID))
    temp=fullLabel(fullLabel(:,2)==0 & fullLabel(:,3)==0 & fullLabel(:,1)==i,6);
    for ii=1:11
        por(ii,1)=sum(temp(:,1)==ii)/size(temp,1);
    end
    wtfcollect(i,:)=por';
end



%% Code Graveyard

% disp([num2str(length(irr(:,2))),' Labeled'])
% disp([num2str(sum(irr(:,2))),' Yes'])
% disp([num2str(length(irr(irr(:,2)==0,1))),' No'])
% disp(['Score: ',num2str(sum(irr(:,2))/length(irr(:,2)))])
% disp(['Most Incorrectly Classified: ',strategylist{mode(irr(irr(:,2)==0,1))}])

% function KEEP(varargin)
% if isempty(varargin);return;end;wh = evalin('caller','who');
% if isempty(wh);error('  There is nothing to keep!');end
% variable = [];for i = 1:length(wh);variable = [variable,':',wh{i}];end;variable = [variable,':'];flag = 0;
% for i = 1:length(varargin)
%     I = findstr(variable,[':',varargin{i},':']);
%     if isempty(I)
%         flag = 1;
%     elseif I == 1
%         variable = variable(1+length(varargin{i})+1:length(variable));
%     elseif I+length(varargin{i})+1 == length(variable)
%         variable = variable(1:I);
%     else
%         variable = [variable(1:I),variable(I+length(varargin{i})+2:length(variable))];
%     end
% end
% if flag == 1;return;end;I = findstr(variable,':');
% if length(I) ~= 1
%     for i = 1:length(I)-1
%         if i ~= length(I)-1;del(i) = {[variable(I(i)+1:I(i+1)-1),' ']};else;del(i) = {variable(I(i)+1:length(variable)-1)};end
%     end
%     evalin('caller',['clear ',del{:}])
% end
% end
