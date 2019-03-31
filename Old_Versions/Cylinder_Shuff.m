% Cylinder_Shuff
%
% Shuffles spike locations and recalculates info content to create a
% shuffled distribution which a 95th percentile can be taken

% Ryan Harvey 9/17

tic
addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analyses\spikeCode'));
addpath(genpath('d:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analyses\spikeCode'));

rats={'LS17','LS19','LS21','LS23'}; %,'LS17','LS19','LS21','LS23'};
datafolder='F:\Users\reharvey\Place_Cell_Data\PAE_Rat\';

rSHUFF=[];
for irats=1:length(rats)
    parent=strcat(datafolder,rats(irats));
    disp(['CYCLING THROUGH RAT:',char(rats(irats))])
    parent=char(parent);
    structdir=dir(parent);
    for I=1:length(structdir) % 1 TO # OF FILES IN DIR
        if structdir(I).isdir && structdir(I).name(1) ~= '.' % IF DIR NOT '.'
            if exist(([parent filesep structdir(I).name filesep 'TT']),'file'); % LOCATE TT FOLDER
                cd([parent filesep structdir(I).name filesep 'TT']); % CD TO TT FOLDER
                CurentFolder=pwd;
                CurrentMat=dir('*.mat'); % LOCATE .MAT FILES IN TT FOLDER
                for J=1:length(CurrentMat) % 1 TO # OF .MAT FILES IN TT FOLDER
                    CurrentMatworking=CurrentMat(J).name;
                    if sum(ismember(CurrentMatworking,'spikeData'))==11 && sum(ismember(CurrentMatworking,'pathProperties'))~=16
                        warning('off','MATLAB:load:variableNotFound');
                        load(CurrentMatworking,'linear_track'); if isequal(linear_track,'yes'); continue; end
                        load(CurrentMatworking,'PeakRate','OverallFR','nSpikes','Field2Wall','data_video_smoothfilt','data_video_smoothfilt3','event');
                        if PeakRate>.25 && nSpikes>50 && OverallFR<=10 && event==2 && exist('Field2Wall','var')
                            disp(['READING: ',CurentFolder filesep CurrentMatworking])
                            spikes=data_video_smoothfilt(:,6);
                            for x = 1:400
                                data_video_smoothfilt(:,6) = circshift(spikes,randi([-size(data_video_smoothfilt,1) size(data_video_smoothfilt,1)],1));
                                [ SmoothRateMap,nBinsx,nBinsy,occ,~] = bindata(data_video_smoothfilt3,30,data_video_smoothfilt(data_video_smoothfilt(:,6) == 1,:),'no',76.5);
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
    keep('rats','irats','parent','datafolder','rSHUFF')
end
NintyFith_Percentile = prctile(rSHUFF,95)

toc



