%BClarkLab Startup Script LB May 2019

%Navigate to folder containing ephys_tools
s=what('ephys_tools');
cd(s.path)

%Initialize Paths 
addpath(fullfile('Analysis'),...
    genpath(fullfile('io')),...
    fullfile('LFP'),...
    genpath(fullfile('preprocessing')),...
    fullfile('Utils'),...
    fullfile('Visualize'),...
    fullfile('external_packages','CircStat2012a'),...
    genpath(fullfile('external_packages','buzcode','externalPackages','FMAToolbox')),...
    fullfile('external_packages','Pass_Index'),...
    fullfile('external_packages','plotSpikeRaster'),...
    fullfile('external_packages','panel'),...
    genpath(fullfile('external_packages','MClust-4.4')),...
    fullfile('external_packages','Colormaps'))

