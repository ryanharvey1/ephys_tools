% postprocessFigures
% plots figures for every cell from standard bclarklab data structure with
% a call to of ( postprocessFigures.main(data) )
%
% You can also call any other of the listed methods shown below independently
%
% List of Methods
%   main
%   raster
%   plot_corr_line
%   unlinearpath
%   spikesonpath
%   spikesonpath_2d
%   ratemaps
%   ratemaps_2d
%   phase_by_pos
%   phase_by_pos_2d
%   phase_map
%   autocors
%   avg_waveforms
%   phase_colormap
%   plot_HD_tuning
%

% Laura Berkowitz March 2019

classdef analyzeLFP_TESTING
    
    methods(Static)
        function p=main(data)
            com=which('postprocessFigures');
            com=strsplit(com,filesep);
            
            basedir=[com{1},filesep,'Users',filesep,com{3},filesep,'GoogleDrive',filesep,'MatlabDir'];
            addpath([basedir,filesep,'BClarkToolbox',filesep,'Analysis',filesep,'Visualize',filesep,'panel'],...
                [basedir,filesep,'BClarkToolbox',filesep,'Analysis'],...
                [basedir,filesep,'BClarkToolbox',filesep,'Analysis',filesep,'Utils'],...
                [basedir,filesep,'CircStat2012a'],...
                [basedir,filesep,'plotSpikeRaster']);
            
            % how many cells, how many sessions?
            nsessions=size(data.events,2);

        end
     
        %function signal_filtered = BandpassFilter(signal, Fs, Fpass)
        
     
    end
end