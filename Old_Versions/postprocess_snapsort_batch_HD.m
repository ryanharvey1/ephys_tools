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

function postprocess_snapsort_batch_HD(Laura,Ryan,figures,rats)
addpath(genpath('D:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis'));
addpath(genpath('d:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis'));

if Ryan==1
    %     rats={'RH11','RH13','RH14','RH16','LS17','LS19','LS21','LS23'};
    datafolder='D:\Place_Cell_Data\RawPAE_PlaceCell\';
    linear_track = 'yes';
    load('D:\Place_Cell_Data\RawPAE_PlaceCell\data.mat','data')
    ratID=fieldnames(data);
    sessions=[];
    for i=1:length(ratID)
        sessions=[sessions;fieldnames(data.(ratID{i}))];
    end
    
elseif Laura==1
    %     rats={'LB01','LB02','LB03','LB04','LB05','LB06'};
    datafolder='F:\P30_Recordings\';
    linear_track = 'no';
    load('F:\P30_Recordings\data.mat','data')
    ratID=fieldnames(data);
    sessions=[];
    for i=1:length(ratID)
        sessions=[sessions;fieldnames(data.(ratID{i}))];
    end
end
h = waitbar(0,'Initializing waitbar...');

% CYCLE THROUGH ALL DATA
for irats=1:length(rats)
    parent = strcat(datafolder,rats(irats)); % CHANGE BASED ON PARENT FOLDER
    disp(['CYCLING THROUGH RAT:',char(rats(irats))])
    parent=char(parent);
    structdir=dir(parent);
    for I=1:length(structdir) % 1 TO # OF FILES IN DIR
        if structdir(I).isdir && structdir(I).name(1) ~= '.'  % IF NOT '.'
            path=[parent filesep structdir(I).name]; % SET PATH TO DATA THAT HAS BEEN THROUGH MCLUST
            cd 'D:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis'
            if sum(ismember(sessions,['S',strjoin(regexp(structdir(I).name,'\d*','Match'),'')]))<1 && exist([path,filesep,'SNAPSorterResults'],'dir')==7
                track_length=TrackLength(path); % SET TRACK LENGTH
                postprocess_snapsort_HDshortVer(path,track_length,linear_track); % CALL POSTPROCESS FUNCTION
            end
        end
    end
%         perc=(irats/length(rats))*100;
%         waitbar(perc/100,h,sprintf('%d%% done...',round(perc)))
end
close(h)
end




