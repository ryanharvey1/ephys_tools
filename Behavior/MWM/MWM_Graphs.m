% MWM_Graphs
close all;clear;clc
load('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/ClarkP30_WaterMaze/resultslbrh8272017.mat')
data=resultslbrh8272017;

% GROUPING INDEX
% WT
WTMale=17:22; WTFemale=23:28;
% TG
TGMale=[1:7,15]; TGFemale=[8:14,16];


% INDEX TRIAL
trial=zeros(length(data),1);
for i=1:length(data)
    trial(i,1)= str2double(regexpi(data{i,1}, '(?<=t\s*)\d*', 'match'));
end

% delete string ID
data=cell2mat(data(:,2:end));

% SEPERATE TRIAL BY GROUP
for i=1:length(unique(trial))
    workingdata=data(trial==i,:);
    WT.(['trial_',num2str(i)])=workingdata(workingdata(:,2)==1,:);
    Tg.(['trial_',num2str(i)])=workingdata(workingdata(:,2)==2,:);
end

% WT Plot
StrategyFig=figure; StrategyFig.Color=[1 1 1];
subplot(1,2,1)
for i=1:length(unique(trial))
    temp=WT.(['trial_',num2str(i)]);
    for ii=1:9
        por(ii,1)=sum(temp(:,3)==ii-1)/size(temp,1);
    end
    collect(i,:)=por';
end
b1=bar(collect,'stacked');
box off
ylim([0 1]); xlim([.5 21.5]);
xlabel('Trial')
ylabel('Proportion of Strategy')
title('WT')
ax1=gca;set(ax1,'FontSize',20,'FontWeight','bold')

%  Tg Plot
subplot(1,2,2)
for i=1:length(unique(trial))
    temp=Tg.(['trial_',num2str(i)]);
    for ii=1:9
        por(ii,1)=sum(temp(:,3)==ii-1)/size(temp,1);
    end
    collect(i,:)=por';
end
b2=bar(collect,'stacked');
box off
ylim([0 1]);xlim([.5 21.5]);
xlabel('Trial')
title('Tg')
ax2=gca;set(ax2,'FontSize',20,'FontWeight','bold')

legend(b2,{'Undefined','Thigmotaxis','Incursion','Scanning','FocusedSearch','ChainingResponse','SelfOrienting','ScanningSurroundings','TargetScanning'},'FontSize',20,'Location','best')

print(StrategyFig, '-dpdf', '-r600','StrategyFig.pdf')





% PLOT GROUP BY SEX~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% SEPERATE TRIAL BY GROUP
for i=1:length(unique(trial))
    workingdata=data(trial==i,:);
    GBS.WTM.(['trial_',num2str(i)])=workingdata(ismember(workingdata(:,1),WTMale),:);
    GBS.TGM.(['trial_',num2str(i)])=workingdata(ismember(workingdata(:,1),TGMale),:);
    GBS.WTF.(['trial_',num2str(i)])=workingdata(ismember(workingdata(:,1),WTFemale),:);
    GBS.TGF.(['trial_',num2str(i)])=workingdata(ismember(workingdata(:,1),TGFemale),:);
end

StrategyGBSFig=figure; StrategyGBSFig.Color=[1 1 1];
subplot(2,2,1)
for i=1:length(unique(trial))
    temp=GBS.WTM.(['trial_',num2str(i)]);
    for ii=1:9
        por(ii,1)=sum(temp(:,3)==ii-1)/size(temp,1);
    end
    collect(i,:)=por';
end
b1=bar(collect,'stacked');
box off
ylim([0 1]); xlim([.5 21.5]);
xlabel('Trial')
ylabel('Proportion of Strategy')
title('WTM')
ax1=gca;set(ax1,'FontSize',20,'FontWeight','bold')

%  Tg Plot
subplot(2,2,2)
for i=1:length(unique(trial))
    temp=GBS.WTF.(['trial_',num2str(i)]);
    for ii=1:9
        por(ii,1)=sum(temp(:,3)==ii-1)/size(temp,1);
    end
    collect(i,:)=por';
end
b2=bar(collect,'stacked');
box off
ylim([0 1]);xlim([.5 21.5]);
xlabel('Trial')
title('WTF')
ax2=gca;set(ax2,'FontSize',20,'FontWeight','bold')

subplot(2,2,3)
for i=1:length(unique(trial))
    temp=GBS.TGM.(['trial_',num2str(i)]);
    for ii=1:9
        por(ii,1)=sum(temp(:,3)==ii-1)/size(temp,1);
    end
    collect(i,:)=por';
end
b1=bar(collect,'stacked');
box off
ylim([0 1]); xlim([.5 21.5]);
xlabel('Trial')
title('TGM')
ax1=gca;set(ax1,'FontSize',20,'FontWeight','bold')

%  Tg Plot
subplot(2,2,4)
for i=1:length(unique(trial))
    temp=GBS.TGF.(['trial_',num2str(i)]);
    for ii=1:9
        por(ii,1)=sum(temp(:,3)==ii-1)/size(temp,1);
    end
    collect(i,:)=por';
end
b2=bar(collect,'stacked');
box off
ylim([0 1]);xlim([.5 21.5]);
xlabel('Trial')
title('TGF')
ax2=gca;set(ax2,'FontSize',20,'FontWeight','bold')

legend(b2,{'Undefined','Thigmotaxis','Incursion','Scanning','FocusedSearch','ChainingResponse','SelfOrienting','ScanningSurroundings','TargetScanning'},'FontSize',20,'Location','best')

% MALE WT
% MALE TG
% FEMALE WT
% FEMALE TG

