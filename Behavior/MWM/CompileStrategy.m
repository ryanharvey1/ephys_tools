%% CompileStrategy
close all
ID=fieldnames(pathResults);
all_strategies=[];
strategyID=[];
loopStats=[];
loopID=[];
loopChi=[];
segID=[];
segmentStats=[];
num_loops=[];
SEX=repmat([1 1 1 1 1 1 1 0 0 0 0 0 0 0 1 0 1 1 1 1 1 0 0 0 0 0 0]',20,1);
GENOTYPE=repmat([1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0]',20,1);
Tg={'sbj1' 'sbj2' 'sbj3' 'sbj4' 'sbj5' 'sbj6' 'sbj7' 'sbj15' 'sbj8' 'sbj9' 'sbj10' 'sbj11' 'sbj12' 'sbj13' 'sbj14' 'sbj16'};
WT={'sbj17' 'sbj18' 'sbj19' 'sbj20' 'sbj21' 'sbj22' 'sbj23' 'sbj24' 'sbj25' 'sbj26' 'sbj27' 'sbj28'};
TgM={'sbj1' 'sbj2' 'sbj3' 'sbj4' 'sbj5' 'sbj6' 'sbj7' 'sbj15'};
TgF={'sbj8' 'sbj9' 'sbj10' 'sbj11' 'sbj12' 'sbj13' 'sbj14' 'sbj16'};
WTM={'sbj17' 'sbj18' 'sbj19' 'sbj20' 'sbj21' 'sbj22'};
WTF={'sbj23' 'sbj24' 'sbj25' 'sbj26' 'sbj27' 'sbj28'};
Male={'sbj1' 'sbj2' 'sbj3' 'sbj4' 'sbj5' 'sbj6' 'sbj7' 'sbj15' 'sbj17' 'sbj18' 'sbj19' 'sbj20' 'sbj21' 'sbj22'};
Female={ 'sbj8' 'sbj9' 'sbj10' 'sbj11' 'sbj12' 'sbj13' 'sbj14' 'sbj16' 'sbj23' 'sbj24' 'sbj25' 'sbj26' 'sbj27' 'sbj28'};
day1={'t1' 't2' 't3' 't4'};
day2={'t5' 't6' 't7' 't8'};
day3={'t9' 't10' 't11' 't12'};
day4={ 't13' 't14' 't15' 't16'};
day5={ 't17' 't18' 't19' 't20'};
NumofStrat=0:11;
for ii=1:length(ID)
    all_strategies=[all_strategies;pathResults.(ID{ii}).strategy];
    strategyID=[strategyID;repmat(ID(ii),length(pathResults.(ID{ii}).strategy),1)];
    
    if sum(isnan(pathResults.(ID{ii}).loopStats))>1
        num_loop=0;
    else
        num_loop=size(pathResults.(ID{ii}).loopStats,1);
    end
    
    
    num_loops=[num_loops;num_loop];
    %     num_loopID=[num_loopID;repmat(ID(ii),length(pathResults.(ID{ii}).strategy),1)];
    loopStats=[loopStats;pathResults.(ID{ii}).loopStats];
    
    
    loopID=[loopID;repmat(ID(ii),size(pathResults.(ID{ii}).loopStats,1),1)];
    
    trial=(strsplit(ID{i,:},'_'));
    strategytemp=pathResults.(ID{ii}).strategy;
    for i=1:length(strategytemp)
        if ismember(trial{3},WT)==1 
           figure(i+1);plot(pathResults.(ID{ii}).segments(i).segments(:,1),pathResults.(ID{ii}).segments(i).segments(:,2)); hold on;
        end
    end

    trial=(strsplit(ID{i,:},'_'));
    strategytemp=pathResults.(ID{ii}).strategy;
    for i=1:length(strategytemp)
           figure(strategytemp(i)+1);plot(pathResults.(ID{ii}).segments(i).segments(:,1),pathResults.(ID{ii}).segments(i).segments(:,2)); hold on;
    end
    
%     for ii=1:length(ID)
%         if pathResults.(ID{ii}).WholeTrjStat(1,5)>= 200
%             for s=1:length(pathResults.(ID{ii}).segments)
%                 if pathResults.(ID{ii}).strategy(s,1)==9 %pathResults.(ID{ii}).segmentStats(s,12)<1 %&& pathResults.(ID{ii}).segmentStats(s,3)<.6 %pathResults.(ID{ii}).strategy(s,1)==6 pathResults.(ID{ii}).strategy(s,1)==6 &&
%                 x=pathResults.(ID{ii}).segments(s).segments(:,1);x=x(:,1);
%                 y=pathResults.(ID{ii}).segments(s).segments(:,2);y=y(:,1);
%                 figure(1); plot(x,y,'LineWidth',1); hold on;
%                 pause(.25)
%                 end
%             end
%         end
%     end
    
%     loopGroup1=[];
%     for i=1:length(ID)
%         trial=(strsplit(ID{i,:},'_'));
%         if ismember(trial{3},TgM); %TGM
%             loopGroup=1;
%             
%         elseif ismember(trial{3},TgF); %TGF
%             loopGroup=2;
%             
%         elseif ismember(trial{3},WTM); %WTM
%             loopGroup=3;
%             
%         else ismember(trial{3},WTF); %WTF
%             loopGroup=4;
%             
%         end
%        
%         loopGroup1=[loopGroup1; loopGroup];
%     end
    segmentStats=[segmentStats;pathResults.(ID{ii}).segmentStats];
    segID=[segID;repmat(ID(ii),size(pathResults.(ID{ii}).segmentStats,1),1)];
end
%% INDEX FOR CHI SQUARE ANALYSIS ON LOOPING
%     for k=1:length(loopID)
%           trial=(strsplit(loopID{k,:},'_'));
%           %genotype
%          if ismember(trial{3},Tg); 
%             loopGChi(k,1)=1;
%          elseif ismember(trial{3},WT); 
%             loopGChi(k,1)=0;
%          end
%          %sex
%         if ismember(trial{3},Male); 
%             loopGChi(k,2)=1;
%          elseif ismember(trial{3},Female); 
%             loopGChi(k,2)=0;
%         end
% %         
%         if ismember(trial{3},TgM); 
%             loopGChi(k,6)=1;
%          elseif ismember(trial{3},TgF); 
%             loopGChi(k,6)=2;
%         elseif ismember(trial{3},WTM); 
%             loopGChi(k,6)=3;
%         elseif ismember(trial{3},WTF); 
%             loopGChi(k,6)=4;
%         end
%         
%         %loop
%         if isnan(loopStats(k,1))==1;
%             loopGChi(k,3)=0;
%         elseif isnan(loopStats(k,1))==0;
%             loopGChi(k,3)=1;
%         end
%         %trial
%         if ismember(trial{2},day1)
%             loopGChi(k,4)=1;
%         elseif ismember(trial{2},day2)
%             loopGChi(k,4)=2;
%         elseif ismember(trial{2},day3)
%             loopGChi(k,4)=3;
%         elseif ismember(trial{2},day4)
%             loopGChi(k,4)=4;
%         elseif ismember(trial{2},day5)
%             loopGChi(k,4)=5;
%         end
%         %trial
%         if strcmp(trial{2},'t1')
%             loopGChi(k,5)=1;
%         elseif strcmp(trial{2},'t2')
%             loopGChi(k,5)=2;
%         elseif strcmp(trial{2},'t3')
%             loopGChi(k,5)=3;
%         elseif strcmp(trial{2},'t4')
%             loopGChi(k,5)=4;
%         elseif strcmp(trial{2},'t5')
%             loopGChi(k,5)=5;
%         elseif strcmp(trial{2},'t6')
%             loopGChi(k,5)=6;
%         elseif strcmp(trial{2},'t7')
%             loopGChi(k,5)=7;
%         elseif strcmp(trial{2},'t8')
%             loopGChi(k,5)=8;
%         elseif strcmp(trial{2},'t9')
%             loopGChi(k,5)=9;
%         elseif strcmp(trial{2},'t10')
%             loopGChi(k,5)=10;
%         elseif strcmp(trial{2},'t11')
%             loopGChi(k,5)=11;
%         elseif strcmp(trial{2},'t12')
%             loopGChi(k,5)=12;
%         elseif strcmp(trial{2},'t13')
%             loopGChi(k,5)=13;
%         elseif strcmp(trial{2},'t14')
%             loopGChi(k,5)=14;
%         elseif strcmp(trial{2},'t15')
%             loopGChi(k,5)=15;
%         elseif strcmp(trial{2},'t16')
%             loopGChi(k,5)=16;
%         elseif strcmp(trial{2},'t17')
%             loopGChi(k,5)=17;
%         elseif strcmp(trial{2},'t18')
%             loopGChi(k,5)=18;
%         elseif strcmp(trial{2},'t19')
%             loopGChi(k,5)=19;
%         elseif strcmp(trial{2},'t20')
%             loopGChi(k,5)=20;
%         end
%     end
    
%     for i=1:27:540
%     loopData3d(:,:,ii)=loopData(i:i+26,:);
%     ii=ii+1;
%     end
%% GROUPING INDEX
% WT
WTMale=17:22; WTFemale=23:28;
% TG
TGMale=[1:7,15]; TGFemale=[8:14,16];
%%
% INDEX TRIAL
trial=zeros(length(strategyID),1);
for i=1:length(strategyID)
    T=str2double(regexpi(strategyID{i,1}, '(?<=t\s*)\d*', 'match'));
    trial(i,:)=T(2); 
end
%%
% SEPERATE TRIAL BY GROUP
group=zeros(length(strategyID),1);
for i=1:length(strategyID)
    T=str2double(regexpi(strategyID{i,1}, '(?<=j\s*)\d*', 'match'));
    group(i,:)=T(1); 
end
group=[group,zeros(length(group),1)];
group(ismember(group(:,1),TGMale),2)=1;
group(ismember(group(:,1),WTMale),2)=2;
group(ismember(group(:,1),TGFemale),2)=3;
group(ismember(group(:,1),WTFemale),2)=4;

strategiestemp=[all_strategies,group(:,2)];

for i=1:length(unique(trial))
    workingdata=strategiestemp(trial==i,:);
    tgm.(['trial_',num2str(i)])=workingdata(workingdata(:,2)==1,:);
    wtm.(['trial_',num2str(i)])=workingdata(workingdata(:,2)==2,:);
    tgf.(['trial_',num2str(i)])=workingdata(workingdata(:,2)==3,:);
    wtf.(['trial_',num2str(i)])=workingdata(workingdata(:,2)==4,:);
end
%%
% WT Plot
StrategyFig=figure; StrategyFig.Color=[1 1 1];
subplot(2,2,1)
for i=1:length(unique(trial))
    temp=tgm.(['trial_',num2str(i)]);
    for ii=1:13
        por(ii,1)=sum(temp(:,1)==ii-1)/size(temp,1);
    end
    tgmcollect(i,:)=por';
end
b1=bar(tgmcollect,'stacked');
box off
ylim([0 1]); xlim([.5 21.5]);
xlabel('Trial')
ylabel('Proportion of Strategy')
title('WT')
ax1=gca;set(ax1,'FontSize',20,'FontWeight','bold')
%
%  Tg Plot
subplot(2,2,2)
for i=1:length(unique(trial))
    temp=wtm.(['trial_',num2str(i)]);
    for ii=1:13
        por(ii,1)=sum(temp(:,1)==ii-1)/size(temp,1);
    end
    wtmcollect(i,:)=por';
end
b2=bar(wtmcollect,'stacked');
box off
ylim([0 1]);xlim([.5 21.5]);
xlabel('Trial')
title('Tg')
ax2=gca;set(ax2,'FontSize',20,'FontWeight','bold')
%  Tg Plot
subplot(2,2,3)
for i=1:length(unique(trial))
    temp=tgf.(['trial_',num2str(i)]);
    for ii=1:13
        por(ii,1)=sum(temp(:,1)==ii-1)/size(temp,1);
    end
    tgfcollect(i,:)=por';
end
b2=bar(tgfcollect,'stacked');
box off
ylim([0 1]);xlim([.5 21.5]);
xlabel('Trial')
title('Tg')
ax2=gca;set(ax2,'FontSize',20,'FontWeight','bold')
%  Tg Plot
subplot(2,2,4)
for i=1:length(unique(trial))
    temp=wtf.(['trial_',num2str(i)]);
    for ii=1:13
        por(ii,1)=sum(temp(:,1)==ii-1)/size(temp,1);
    end
    wtfcollect(i,:)=por';
end
b2=bar(wtfcollect,'stacked');
box off
ylim([0 1]);xlim([.5 21.5]);
xlabel('Trial')
title('Tg')
ax2=gca;set(ax2,'FontSize',20,'FontWeight','bold')

legend(b2,{'Undefined','directswim','directed swim','Circuitous Direct','TargetSearch','TargetScanning','Chaining','focused search','scanning','scanning surround','incursion','self orienting','thigmotaxis'},'FontSize',20,'Location','best')

% GROUPING INDEX
% WT
WTMale=17:22; WTFemale=23:28;
% TG
TGMale=[1:7,15]; TGFemale=[8:14,16];
%%
% INDEX TRIAL
trial=zeros(length(strategyID),1);
for i=1:length(strategyID)
    T=str2double(regexpi(strategyID{i,1}, '(?<=t\s*)\d*', 'match'));
    trial(i,:)=T(2); 
end
%%
% SEPERATE TRIAL BY GROUP
group=zeros(length(strategyID),1);
for i=1:length(strategyID)
    T=str2double(regexpi(strategyID{i,1}, '(?<=j\s*)\d*', 'match'));
    group(i,:)=T(1); 
end
group=[group,zeros(length(group),1)];
group(ismember(group(:,1),[WTMale,WTFemale]),2)=1;
group(ismember(group(:,1),[TGMale,TGFemale]),2)=2;

strategiestemp=[all_strategies,group(:,2)];

for i=1:length(unique(trial))
    workingdata=strategiestemp(trial==i,:);
    wt.(['trial_',num2str(i)])=workingdata(workingdata(:,2)==1,:);
    tg.(['trial_',num2str(i)])=workingdata(workingdata(:,2)==2,:);
end
%%
% WT Plot
StrategyFig=figure; StrategyFig.Color=[1 1 1];
subplot(1,2,1)
for i=1:length(unique(trial))
    temp=wt.(['trial_',num2str(i)]);
    for ii=1:15
        por(ii,1)=sum(temp(:,1)==ii-1)/size(temp,1);
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
%
%  Tg Plot
subplot(1,2,2)
for i=1:length(unique(trial))
    temp=tg.(['trial_',num2str(i)]);
    for ii=1:15
        por(ii,1)=sum(temp(:,1)==ii-1)/size(temp,1);
    end
    collect(i,:)=por';
end
b2=bar(collect,'stacked');
box off
ylim([0 1]);xlim([.5 21.5]);
xlabel('Trial')
title('Tg')
% ax2=gca;set(ax2,[.7 .7 .7],'FontSize',20,'FontWeight','bold')

legend(b2,{'Undefined','directswim','selforienteddirectedswim','directed swim','Arching Direct','Circuitous Direct','Thigmotaxis','scanning','FocalSearch','TargetSearch','OppTargetSearch','SelfOrienting','Chaining','Incursion','scanningsurrounding','targetscanning'},'FontSize',20,'Location','best')

% print(StrategyFig, '-dpdf', '-r600','StrategyFig.pdf')


%%



% MALE WT
% MALE TG
% FEMALE WT
% FEMALE TG







% % histogram(all_strategies)
% 
% % PLOT ALL STRATEGIES
% NumofStrat=unique(all_strategies);
% for ii=0:length(NumofStrat)
%     [I]=find(all_strategies==ii);
%     for i=1:length(I)
%         figure(ii+1)
%         plot(classification_configs.CLASSIFICATION.segments.items(1, I(i)).points(:,2),classification_configs.CLASSIFICATION.segments.items(1, I(i)).points(:,3))
%         hold on
%     end
% end
% 
% 
% 
% all_strategiesSort=sort(all_strategies)