% ClustSum1Script  Script to save Matlab figure files for "Cluster Summary" (BPALL) sheets
%
%  - must have .tt (or .dat if data were upsampled) files in appropriate 'TT' folders to work
%  - modify this script to be specific to whatever 'TT' folders exist for each session, and to whatever
%     the .t and .dat files are named in a given folder
%  - start with Matlab in the first 'TT' folder in the list below
%
% INPUTS to check_cluster_from_file2 commands below:
%   epochs - a struct object with fields:
%               epochs.names = cell array of strings of epoch names 
%               epochs.intervals = cell array of 1x2 arrays with [start_ts  end_ts] (start
%                   and end timestamps of each epoch) -- elements in cell array correspond
%                   to epochs, listed in same sequence as in epochs.names
%   2nd input ('TT2', etc.) - .t file naming convention for given tetrode (e.g., 'TTctt_' or 'TT_')
%   3rd input ('tt2.dat', etc.) - name of .tt file, or .dat file if upsampled, for given tetrode
%   4th & 5th inputs - use 'save_figures' and 'fig' to save Matlab .fig files (one per
%       cluster) in 'TT' folders
%
% MRN 2/04


if exist('epochs.mat')
    load epochs;
end %if exist('epochs.mat')

TTfolders = dir;
for i=3:length(TTfolders)
    if TTfolders(i).isdir
        if ( strcmp(TTfolders(i).name(1:2),'TT') | strcmp(TTfolders(i).name(1:2),'tt') )

            cd(TTfolders(i).name)
            
            tfiles = findfiles('*.t');
            [d, fn, ext] = fileparts(tfiles{1});
            fname = fn(1:find(fn == '_')-1);
            
            if exist([fname '.tt'])
                tt_or_dat = [fname '.tt'];
            else
                if exist([fname '.dat'])
                    tt_or_dat = [fname '.dat'];
                else 
                    if exist([fname '.DAT'])
                        tt_or_dat = [fname '.DAT'];
                    else
                        if exist([fname '.ntt'])
                            tt_or_dat = [fname '.ntt'];
                        end %if exist([fname '.ntt'])
                    end %if exist([fname '.DAT'])
                end %if exist([fname '.dat'])
            end %if exist([fname '.tt'])
            
            ClustSum1(epochs,fname,tt_or_dat,'save_figures','fig');
            
            cd ..
            
        end %if ( strcmp ...
    end %if TTfolders ...
end %for

clear TTfolders i tfiles d fn ext fname tt_or_dat;