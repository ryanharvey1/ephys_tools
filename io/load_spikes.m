function [S,avgwave,ID,cellnum,tetrode,clusterquality]=load_spikes(path)
% load_spikes: loads SNAPSorter* format spike times and spike info 
%
%   Input
%           path: path to spike data (TT*.ntt, TT*_info.mat, & TT*.mat files
%                   created by SNAPsorter or after_spikesort_cleanup.m)
%   Output
%           S: cell array of spike times in microseconds
%           avgwave: average waveforms 
%           ID: full path to cell
%           cellnum: cell number
%           tetrode: tetrode number
%
% *After testing, I really don't recommend using SNAPSorter
%
% Ryan Harvey

% sets path to snap folder & looks within for mat files
if ~exist([path,filesep,'Sorted'],'dir')
    error('Make sure to run after_spikesort_cleanup.m and cd to session dir')
end
sortedpath = [path,filesep,'Sorted'];
file=struct2table(dir('Sorted/*.mat'));
t=table2cell(file(:,1));
file=t(~contains(t,'_info.mat'),1);
info=strcat(erase(file,'.mat'),repmat('_info.mat',length(file),1));

% EXTRACT SPIKE TIMES,ID, AND AVERAGE WAVEFORMS FROM SNAPSORT .MAT FILES
ns=1;clusterquality=[];
for m=1:length(file)
    load(fullfile(sortedpath,file{m}));
    load(fullfile(sortedpath,info{m}));
    clusterquality=[clusterquality;grades(:,[1,3,5])];
    cells=unique(output(:,2));
    idxmeans=1;
    for n=1:length(cells)
        tetrode{ns,1}=file{m};
        cellnum(ns,1)=cells(n);
        ID{ns,1}=orig_filename;
        avgwave{ns,1}=means{1,idxmeans};
        S{ns,1}=output(output(:,2)==cells(n),1);
        ns=ns+1;
        idxmeans=idxmeans+1;
    end
end

% Snap sorter outputs ts in ms. We need to convert to microseconds
for spikes=1:length(S)
    if sum(S{spikes}<1e8)==length(S{spikes})
        S{spikes}=S{spikes}*1000000;
    end
end
end