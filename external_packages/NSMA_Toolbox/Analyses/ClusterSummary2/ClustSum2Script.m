% ClustSum2Script  Runs the ClustSum2 function for all cells from a given session
%
% Produces and saves "Cluster Summary 2" figures for all t-files in a given session.
% Run from folder containing "TT" subdirectories (.t files within these subdirs).
% Make sure 'epochs.mat' and 'position.mat' saved workspaces are in this folder, or
%   "epochs" and "X","Y" position data variables are in current workspace. 
% This script should also be in this folder, or in a folder in your path.
% The functions 'ClustSum2', 'TuningCurves', 'findfiles', and 'LoadSpikes'
%   should be in your path.
% Variables S (cell array of ts objects of spikes in each cluster) & tfiles (cell
%   array of t-file names) retained in workspace in case you want to save them.
%
% PL,MN 2002, last modified '04 by MN


% get position data
if exist('position.mat')
    load position X Y;
end %if exist('position.mat')


% get epochs
if exist('epochs.mat')
    load epochs;
end %if exist('epochs.mat')


% get spike data
TTfolders = dir;
tfiles = {};
S = {};
for i=3:length(TTfolders)
    if TTfolders(i).isdir
        if ( strcmp(TTfolders(i).name(1:2),'TT') | strcmp(TTfolders(i).name(1:2),'tt') )

            cd(TTfolders(i).name)
            
            tfiles_temp = findfiles('*.t'); 
            S_temp = LoadSpikes(tfiles_temp);

            tfiles = [tfiles; tfiles_temp];
            S = [S; S_temp];
            
            cd ..
            
        end %if ( strcmp ...
    end %if TTfolders ...
end %for


% Produce ClustSum2 figs & save in current folder
for iCell = 1:length(S)
    
    
    % IF YOU WANT ONLY CERTAIN EPOCHS SHOWN IN SCATTERPLOTS OF SPIKE LOCATIONS:
    epochs_plotted = epochs.names(2:end-1); %this uses only the 2nd, 3rd, ..., 2nd-to-last
                                            % epochs only - ignores 1st & last - this
                                            % assumes that the 1st & last are the sleeps.
    fh = ClustSum2(tfiles{iCell}, S{iCell}, X, Y, epochs, iCell, epochs_plotted);
    
    
    % OR, IF YOU WANT DATA FOR THE ENTIRE SESSION SHOWN:
    %fh = ClustSum2(tfiles{iCell}, S{iCell}, X, Y, epochs, iCell);
    
    
    % save figures as Matlab figure files
    [d, fname, ext] = fileparts(tfiles{iCell});
    fnstrg = ['ClustSum2_' fname];
    saveas(fh,[fnstrg '.fig']);
    
    close all;
end


clear TTfolders i tfiles_temp S_temp iCell fh d fname ext fnstrg;


