
%%
%%%%%%%%%%% Plot Paths with home base outline and cue location
% load('params_V8');
load('params_V14');
param_idx=params.PathName;

vars=fieldnames(params);

%% Results Fig. 1 (Path Length, Search Area, Time in outer wall, running speed)
[row,~]=find(ismember(vars,{'pathL','searchArea','dwellOutsideWall','runSpeed'}));

idx=1;
fig=figure; fig.Color=[1 1 1];

for i=1:size(row,1)
tg1=vertcat(params.(vars{row(i)}){contains(param_idx,'Tg') & contains(param_idx,'D1')});
wt1=vertcat(params.(vars{row(i)}){contains(param_idx,'WT') & contains(param_idx,'D1')});
tg2=vertcat(params.(vars{row(i)}){contains(param_idx,'Tg') & contains(param_idx,'D2')});
wt2=vertcat(params.(vars{row(i)}){contains(param_idx,'WT') & contains(param_idx,'D2')});

subplot(size(row,1),2,idx)
plotspread_wrapper(wt1,tg1,{'WT','Tg'})
ylim([0 max([wt2;tg2;wt1;tg1])])
ylabel(vars{row(i)})
title('Proximal Cue Present')

idx=idx+1;
subplot(size(row,1),2,idx)
plotspread_wrapper(wt2,tg2,{'WT','Tg'})
ylim([0 max([wt2;tg2;wt1;tg1])])
title('Proximal Cue Absent')
idx=idx+1;
end

%% Figure 2 

