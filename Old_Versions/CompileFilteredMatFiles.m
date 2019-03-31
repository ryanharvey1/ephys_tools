%% CompileFilteredMatFiles
% This script will cycle through a parent directory and pull out info from the .mat files
% Ryan Harvey 12/16
clear,close all
if ismac==0
    addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analyses\spikeCode'));
    addpath(genpath('d:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analyses\spikeCode'));
    addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis'));
    
    for igroup=1:2
        PopSmoothRateMap=[];
        SmoothRateMap_Rightw=[];
        SmoothRateMap_Leftw=[];
        if igroup==1
            rats={'RH13','RH14'};
        elseif igroup==2
            rats={'RH11','RH16'};
        end
        for irats=1:length(rats)
            parent=strcat('F:\Users\reharvey\Place_Cell_Data\PAE_Rat\',rats(irats));
            disp(['CYCLING THROUGH RAT:',char(rats(irats))])
            parent=char(parent);
            structdir=dir(parent);
            %         K=1;
            for I=1:length(structdir) % 1 TO # OF FILES IN DIR
                if structdir(I).isdir && structdir(I).name(1) ~= '.' % IF DIR NOT '.'
                    if exist(([parent filesep structdir(I).name filesep 'TT']),'file'); % LOCATE TT FOLDER
                        CurentFolder=cd([parent filesep structdir(I).name filesep 'TT']); % CD TO TT FOLDER
                        CurrentMat=dir('*.mat'); % LOCATE .MAT FILES IN TT FOLDER
                        for J=1:length(CurrentMat) % 1 TO # OF .MAT FILES IN TT FOLDER
                            CurrentMatworking=CurrentMat(J).name;
                            if sum(ismember(CurrentMatworking,'spikeData'))==11 && sum(ismember(CurrentMatworking,'pathProperties'))~=16
                                warning('off','MATLAB:load:variableNotFound');
                                load(CurrentMatworking,'SmoothRateMap_Right','SmoothRateMap_Left','SmoothRateMap','InformationContent','PeakRate','OverallFR','nSpikes'); % OPEN .MAT WITH THESE VARS
                                disp(['READING: ',CurentFolder filesep CurrentMatworking])
                                if exist('SmoothRateMap_Right','var')
                                    if InformationContent>=0.40 && nSpikes>=50 && PeakRate>=.5 && OverallFR<=10
                                        
                                        %SmoothRateMap_Rightw,num_spikes_Rightw]=RightVsLeft(right,SpikeFile,track_length,sampleRate);
                                        if size(SmoothRateMap_Right,2)<22; SmoothRateMap_Right=imresize(SmoothRateMap_Right,[1,22]);end
                                        SmoothRateMap_Rightw=[SmoothRateMap_Rightw;SmoothRateMap_Right];
                                        %num_spikes_Right=[num_spikes_Right;num_spikes_Rightw];
                                        
                                        %[SmoothRateMap_Leftw,num_spikes_Leftw]=RightVsLeft(left,SpikeFile,track_length,sampleRate);
                                        if size(SmoothRateMap_Left,2)<22; SmoothRateMap_Left=imresize(SmoothRateMap_Left,[1,22]);end
                                        SmoothRateMap_Leftw=[SmoothRateMap_Leftw;SmoothRateMap_Left];
                                        %num_spikes_Left=[num_spikes_Left;num_spikes_Leftw];
                                        
                                        %                                     Overall_DirectionalityIndex(K,:)=[num_spikes_Leftw num_spikes_Rightw];
                                        
                                        if size(SmoothRateMap,2)<22; SmoothRateMap=imresize(SmoothRateMap,[1,22]);end
                                        PopSmoothRateMap=[PopSmoothRateMap;SmoothRateMap];
                                        %                                     K=K+1;
                                    end
                                end
                                if ~exist('SmoothRateMap_Right','var') % IF THE .MAT FILE WAS INCORRECT
                                    keep('I','J','K','structdir','parent','CurrentMat','CurrentFolder','rats','irats','igroup');
                                else
                                    keep('I','J','K','SmoothRateMap_Rightw','SmoothRateMap_Leftw','Overall_DirectionalityIndex','PopSmoothRateMapAll','structdir','parent','CurrentMat','CurentFolder','rats','irats','igroup'); % KEEP ALLSPIKEDATA & STRUCTDIR IF .MAT WAS CORRECT
                                end
                            end
                        end
                        keep('I','K','SmoothRateMap_Rightw','SmoothRateMap_Leftw','Overall_DirectionalityIndex','PopSmoothRateMapAll','structdir','parent','CurrentMat','rats','irats','igroup'); % DON'T KEEP VAR J BECAUSE YOU JUST EXITED THE J LOOP
                    end
                end
            end
            disp(['DONE WITH:',char(rats(irats))])
            keep('rats','irats','parent','igroup','SmoothRateMap_Rightw','SmoothRateMap_Leftw','Overall_DirectionalityIndex','PopSmoothRateMapAll')
        end
        if igroup==1
            Group=('Control');
        elseif igroup==2
            Group=('PAE');
        end
        % TAKE FILTERED COMPILED DATA AND PLOT THEM FOR EACH ANIMAL IN EACH DIRECTION
        DirectionalMatrix(parent,[],0,SmoothRateMap_Rightw,SmoothRateMap_Leftw,Group);
        % CREATE POPULATION MAPS
        %     DirectionalMatrix(parent,[],0,PopSmoothRateMap,Group);
        keep('igroup','parent')
    end
    % MAKE FIGURES AND COMPUTE STATS
    %
    [filepath, ~] = fileparts(parent);
else
    addpath(genpath('/Users/ryanharvey/GoogleDrive/MatlabDir/NSMA_Toolbox'))
    addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\NSMA_Toolbox'))
    for igroup=1:2
        if igroup==1
            Group=('Control');
        elseif igroup==2
            Group=('PAE');
        end
        if ismac
            load(['/Users/RyanHarvey/GoogleDrive/postprocessMclust_test' filesep Group '_Data.mat'])
        else
            load([filepath filesep Group '_Data.mat'])
        end
        
        % REMOVE NANS
        SmoothRateMap_Right_arrangedR(isnan(SmoothRateMap_Right_arrangedR))=0;
        SmoothRateMap_Left_arrangedR(isnan(SmoothRateMap_Left_arrangedR))=0;
        
        % AUTO CORRELATE & CROSS CORRELATE
        autoR=corrcoef_AB(SmoothRateMap_Right_arrangedR,SmoothRateMap_Right_arrangedR);
        R=corrcoef_AB(SmoothRateMap_Right_arrangedR,SmoothRateMap_Left_arrangedR);
        
        % FIND DECORRELATING DISTANCE
        decordist=zeros(size(R,2),1);
        for a=0:size(R,2)-1
            decordist(a+1)=nanmean(diag(R, a));
        end
        reversedecorr=flipud(decordist); reversedecorr(end)=[];
        decordistALL(:,igroup)=[reversedecorr;decordist];
        
        %PLOT POPULATION VECTORS
        matfordis=[autoR,R;R,autoR];
        figure(5); subplot(1,2,igroup); h = pcolor(matfordis);
        colormap jet
        axis off
        hold on
        colorbar
        box off
        set(h, 'EdgeColor', 'none');
        set(gca,'YDir','reverse');
        title(['Population vector correlation: ',num2str(Group)]);
    end
    % PLOT DECORRELATING DISTANCE
    figure(6);
    h=plot(linspace(-120,120,43),decordistALL(:,1));
    hold on
    plot(linspace(-120,120,43),decordistALL(:,2));
    xlabel('Distance (cm)')
    ylabel('Correlation')
    title('Decorrelation Distance')
    legend({'Control','PAE'},'FontSize',12)
end
% STATS
% BASE STATS OFF OF THESE
% Kjelstrup... Moser 2008; Battaglia, SUTHERLAND, MCNAUGHTON 2004;
% Ravassard... MEHTA 2013; Maurer...McNaughton 2005

