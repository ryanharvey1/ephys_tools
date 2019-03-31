% postprocessMClust_v9_batch
%
% THIS VERSION CYCLES THROUGH ALL DATA
% THIS SCRIPT CREATES FIGURES AND STATS FOR PLACE & HD CELL DATA
% CAN HANDLE CIRCULAR ARENA & LINEAR TRACK DATA
%
% INPUT:
%       postprocessMClust_v9_batch(Laura,Ryan,figures)
%       postprocessMClust_v9_batch(1,0,true) >>> for Laura's data
%       postprocessMClust_v9_batch(0,1,true) >>> for Ryan's data
%
% OUTPUTS: 
%        RATE MAP / SPIKE ON PATH / TUNING CURVE / POLAR PLOT / R VS. L SPIKES / STATS etc*.
%        *Everything in the normal postprocessing script
%
% NOTE: THIS FUNCTION CALLS COMPILE ALL MAT FILES WHEN FINISHED
% Ryan Harvey & Ben Clark 2016-2017

function postprocessMClust_v9_batchTestcopy(Laura,Ryan,figures)
addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis'));
addpath(genpath('d:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis'));
tic
if Ryan==1
    rats={'RH11','RH13','RH14','RH16','LS17','LS19','LS21','LS23'}; 
    datafolder='/Users/ryanharvey/OneDrive - University of New Mexico/Test_Sessions/';
    linear_track = 'yes'; 
elseif Laura==1
    rats={'LB01','LB02','LB03','LB04','LB05','LB06'};
    datafolder='D:\ClarkP30_Recordings\';
    linear_track = 'no'; 
end
poolobj = gcp;
addAttachedFiles(poolobj,{'postprocessMClust_v9.m'})
% CYCLE THROUGH ALL DATA
parfor irats=1:length(rats)
    parent = strcat(datafolder,rats(irats)); % CHANGE BASED ON PARENT FOLDER
    disp(['CYCLING THROUGH RAT:',char(rats(irats))])
    parent=char(parent);
    structdir=dir(parent);
    for I=1:length(structdir) % 1 TO # OF FILES IN DIR
        if structdir(I).isdir && structdir(I).name(1) ~= '.' && any(regexp(structdir(I).name,'p$'))~=1; % IF NOT '.' & NOT P FILE
            if exist([parent filesep structdir(I).name,'p',filesep 'TT'],'dir');
                cd([parent filesep structdir(I).name,'p',filesep 'TT']); % CD TO TT FOLDER
            end
            if any(size(dir([parent filesep structdir(I).name,'p',filesep 'TT' filesep '*.t' ]),1))==1 % CHECK FOR MCLUST FILES
                path=[parent filesep structdir(I).name]; % SET PATH TO DATA THAT HAS BEEN THROUGH MCLUST
            else
                continue
            end
        else
            continue
        end
        track_length=TrackLength(path); % SET TRACK LENGTH       
        postprocessMClust_v9(path,track_length,linear_track,figures); % CALL POSTPROCESS FUNCTION
    end
end
disp('DONE WITH POST PROCESS')
toc
% disp('STARTING TO COMPILE')
% CompileAllMatFilesv5(Laura,Ryan) % CALL COMPILE ALL MATFILES





