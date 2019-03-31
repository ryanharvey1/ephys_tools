%% Load data
tilteddata=xlsread ('Tilted place cell data.xlsx');
controldata=xlsread ('control place cell data.xlsx');
% filter durations less than 5 out
x=find(tilteddata(:,5)>5);
tilteddata=tilteddata(x,:);
% waveform length by average firing rate
% tilted
figure (1)
tilted=scatter(tilteddata(:,5),tilteddata(:,4),'filled','R');
hold on
% control
figure (1)
control=scatter(controldata(:,5),controldata(:,4),'filled','K');
% Adding labels and label sizes
title=title('Waveform lengths by Avg Firing Rate');
set(title, 'FontSize', 14)
xlabel=xlabel('Waveform Length');
set(xlabel, 'FontSize', 14)
ylabel=ylabel('Avg Firing Rate');
set(ylabel, 'FontSize', 14)
legend=legend([tilted,control],{'Tilted', 'Control'},'location','best');
set(legend, 'FontSize', 14)
% ttest
[h,p,ci,stats] = ttest2(tilteddata(:,5),controldata(:,5))
clear all
%% Graphs seperated out by tilted
% Load data
tilteddata=xlsread ('Tilted place cell data.xlsx');
% filter durations less than 5 out
x=find(tilteddata(:,5)>5);
tilteddata=tilteddata(x,:);
% waveform length by average firing rate
% tilted
figure (2)
tilted=scatter(tilteddata(:,5),tilteddata(:,4),'filled','K');
% Adding labels and label sizes
title=title('Waveform lengths by Avg Firing Rate Tilted');
set(title, 'FontSize', 14)
xlabel=xlabel('Waveform Length');
set(xlabel, 'FontSize', 14)
ylabel=ylabel('Avg Firing Rate');
set(ylabel, 'FontSize', 14)
clear all
% Graphs seperated out by control
controldata=xlsread ('control place cell data.xlsx');
figure (3)
control=scatter(controldata(:,5),controldata(:,4),'filled','K');
% Adding labels and label sizes
title=title('Waveform lengths by Avg Firing Rate Control');
set(title, 'FontSize', 14)
xlabel=xlabel('Waveform Length');
set(xlabel, 'FontSize', 14)
ylabel=ylabel('Avg Firing Rate');
set(ylabel, 'FontSize', 14)
clear all