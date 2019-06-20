
%%
%%%%%%%%%%% Plot Paths with home base outline and cue location
load('params_V14');
param_idx=params.PathName;

vars=fieldnames(params);

%% Results section 1: Gross Measures of Exploratory analysis 

[row,~]=find(ismember(vars,{'pathL','searchArea','dwellOutsideWall','runSpeed'}));

fig=figure; fig.Color=[1 1 1];

for i=1:size(row,1)
tg1=vertcat(params.(vars{row(i)}){contains(param_idx,'Tg') & contains(param_idx,'D1')});
wt1=vertcat(params.(vars{row(i)}){contains(param_idx,'WT') & contains(param_idx,'D1')});
tg2=vertcat(params.(vars{row(i)}){contains(param_idx,'Tg') & contains(param_idx,'D2')});
wt2=vertcat(params.(vars{row(i)}){contains(param_idx,'WT') & contains(param_idx,'D2')});

RL_anova([tg1 tg2],[wt1 wt2],vars{row(i)})

%means
disp(['Tg Day 1: ',vars{i},' mean=',num2str(nanmean(tg1)),' SEM=',num2str(nanstd(tg1)/sqrt(size(tg1,1)))]); 
disp(['WT Day 1: ',vars{i},' mean=',num2str(nanmean(wt1)),' SEM=',num2str(nanstd(wt1)/sqrt(size(wt1,1)))]); 
disp(['Tg Day 2: ',vars{i},' mean=',num2str(nanmean(tg2)),' SEM=',num2str(nanstd(tg2)/sqrt(size(tg2,1)))]); 
disp(['Wt Day 2: ',vars{i},' mean=',num2str(nanmean(wt2)),' SEM=',num2str(nanstd(wt2)/sqrt(size(wt2,1)))]);  

disp(['Overall Tg : ',vars{i},' mean=',num2str(nanmean([tg1;tg2])),' SEM=',num2str(nanstd(tg1)/sqrt(size([tg1;tg2],1)))]); 
disp(['Overall WT : ',vars{i},' mean=',num2str(nanmean([wt1;wt2])),' SEM=',num2str(nanstd(wt1)/sqrt(size([wt1;wt2],1)))]); 


end


