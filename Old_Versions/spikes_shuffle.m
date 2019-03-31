% spikes_shuffle shifts the spike series relative to the time series
%
% Input:
%
% Output:
%       - rSHUFF = Spatial Information Content for shuffled spike series
%
% Created by Ben C 2014; Modified Ben C Feb 03 2017 ; Modified Ryan H 2/5/17

% identify path to data files

clc, clear, close all

addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analyses\spikeCode'));
addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\CircStat2012a'));
addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis'));

rSHUFF = [];
rats={'RH11','RH13','RH14','RH16','LE2813','LE2821'};

tic;
for irats=1:length(rats)
    parent=strcat('F:\Users\reharvey\Place_Cell_Data\PAE_Rat\',rats(irats));
    disp(['CYCLING THROUGH RAT:',char(rats(irats))])
    parent=char(parent);
    structdir=dir(parent);
    K=1;
    for I=1:length(structdir) % 1 TO # OF FILES IN DIR
        if structdir(I).isdir && structdir(I).name(1) ~= '.' % IF DIR NOT '.'
            if exist(([parent filesep structdir(I).name filesep 'TT']),'file'); % LOCATE TT FOLDER
                cd([parent filesep structdir(I).name filesep 'TT']); % CD TO TT FOLDER
                CurrentMat=dir('*.mat'); % LOCATE .MAT FILES IN TT FOLDER
                for J=1:length(CurrentMat) % 1 TO # OF .MAT FILES IN TT FOLDER
                    CurrentMatworking=CurrentMat(J).name;
                    if sum(ismember(CurrentMatworking,'spikeData'))==11 && sum(ismember(CurrentMatworking,'pathProperties'))~=16
                        load(CurrentMatworking,'data_video_smoothfilt','data_video_smoothfilt3','track_length2','PerSpkslessthan2ms','PeakRate')
                        if sum(data_video_smoothfilt(:,6))<50; continue; end % Don't calculate for <50 spikes (Would inflate IC)
                        disp(['READING: ',pwd filesep CurrentMatworking])
                        if exist('data_video_smoothfilt','var')==1 && exist('track_length2','var')==1 && PerSpkslessthan2ms<1.5 && PeakRate>.25 && (sum(data_video_smoothfilt(:,6))/size(data_video_smoothfilt3,1))*30<=10
                            spksf=data_video_smoothfilt(:,6);
                            for x = 1:400;
                                shuff_interv = 30*20;
                                spkSHUFFi = circshift(spksf,randi([shuff_interv size(data_video_smoothfilt,1)-shuff_interv],1));
                                data_video_smoothfilt(:,6)=spkSHUFFi;
                                spks_VEL = data_video_smoothfilt(data_video_smoothfilt(:,6) == 1,:);
                                [ SmoothRateMap,nBinsx,nBinsy,occ,~] = bindata(data_video_smoothfilt3,30,spks_VEL,'yes',track_length2);
                                rY = reshape(SmoothRateMap,nBinsx*nBinsy,1); % reshape data into column
                                rY(isnan(rY))=0; rY(isinf(rY))=0;
                                occRSHP = reshape(occ,nBinsx*nBinsy,1); % reshape data into column
                                occSUM = sum(occRSHP); % summed occupancy
                                pX = occRSHP./occSUM; % normalized occupancy
                                [InformationContent] = InformationPerSpike(rY,pX); % from
                                rSHUFF = [rSHUFF; InformationContent];
                            end
                        end
                    end
                    keep('I','J','K','structdir','parent','CurrentMat','rats','irats','rSHUFF'); % KEEP ALLSPIKEDATA & STRUCTDIR IF .MAT WAS CORRECT
                end
                clear J % DON'T KEEP VAR J BECAUSE YOU JUST EXITED THE J LOOP
            end
        end
    end
    disp(['DONE WITH:',char(rats(irats))])
    keep('rats','irats','parent','rSHUFF')
end
NintyFith_Percentile = prctile(rSHUFF,95)
toc
disp('DONE')



%         clear all non essential variables
%         keep('rSHUFF', 'corrSHUFF', 'path', 'ReadData', 'sampleRate', 'xmin', 'xmax', 'ymin', 'ymax', 'minLED', 'maxLED', 'dBins', 'shuff_interv', 'i');
%         [filepath, filename] = fileparts(ReadData{i});
%         save([filepath filesep filename '_SHUFF.mat']);
% end






