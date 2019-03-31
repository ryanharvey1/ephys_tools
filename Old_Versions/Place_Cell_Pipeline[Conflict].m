%% Place_Cell_Pipeline

% Description: 
% This pipleline was designed to analyze place cell data from either a linear track or an arena
% To go through each step, click anywhere in section of code and click "Run Section" above
% Ryan Harvey & Ben Clark 2016

%% STEP #1      ***PREPROCESS DATA***
% input path to folder or use folder gui to select folder
tic; clear 
addpath(genpath('d:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis'));
addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis'));

% SPECIFY DIRECTOR AND TEMPLATE ('hipp' or 'pfc')              
multible_sessions=true; 
PreprocessingMClust('D:\ClarkP30_Recordings\LB03\2017-04-20_11-50-49','pfc',multible_sessions);

disp 'DONE'; toc % displays time to completion 


%% STEP #2      ***USE MCLUST TO SPIKE SORT DATA***
clear
addpath(genpath('d:\Users\BClarkLab\GoogleDrive\MatlabDir\MClust3.1UofL_c'));
addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\MClust3.1UofL_c'));
%
MClust
    % Basic Instructions:
    % Click "run klustakwik"
    % click "create/load FD..."
    % select NTT file
    % export all clusters in the new window
    % unclick "run klustakwik"
    % click "create/load FD..." again
    % Select NTT file
    % another window will pop up. Click on cluster file you just created
    
    % For help with Cutting/merging with MClust: 
            % Search for the following doc or pdf ('Preparing Data for Cluster Cutting')
            
%% STEP #3      *** POSTPROCESS DATA ***
%
clear; close all
warning('off', 'MATLAB:dispatcher:nameConflict')
addpath(genpath('d:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis'));
addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis'));

tic
postprocessMClust_v9('/Users/RyanHarvey/GoogleDrive/postprocessMclust_test/2016-07-02_18-35-24',120,'yes',true)  
toc
    % 'path' is the path to the non-p folder
    % 'track_length' is the length of the track in centimeters: 90 (pre 7/14/16) or 120 (7/14/16 and on)
    % 'Linear_track' is if the maze was a linear track or arena: if linear track ('yes'), if not ('no')
    
% Example input:  
        % postprocessMClust('F:\Users\reharvey\Place_Cell_Data\PAE_Rat\LS23\2017-03-16_17-29-05','120', 'yes',true)   
%% STEP #4     ***BATCH POSTPROCESSING DATA***

% This step post processes all data
% Only run after major changes to normal post processing function

clear, clc , close all
tic
% CHOOSE WHICH GROUP TO RUN (RYAN OR LAURA)
% Ryan=1;
% Laura=0;
% 
% if Ryan==1
%     rats={'RH11','RH13','RH14','RH16','LS17','LS19','LS21','LS23'}; 
%     datafolder='F:\Users\reharvey\Place_Cell_Data\PAE_Rat\';
%     linear_track = 'yes'; 
% elseif Laura==1
%     rats={'LB01','LB02','LB03','LB04','LB05','LB06'};
%     datafolder='D:\ClarkP30_Recordings\';
%     linear_track = 'no'; 
% end
% parfor ratsi=1:length(rats)
%     postprocessMClust_v9_batch(rats(ratsi),datafolder,linear_track,true)
% end
postprocessMClust_v9_batch(0,1,true,{'RH16'})

toc
% Example
% postprocessMClust_v9_batch(0,1,true)




