% CompileFilteredMatFiles
% This script will cycle through a parent directory and pull out info from the .mat files
% Ryan Harvey 12/16
clear,close all
addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analyses\spikeCode'));
addpath(genpath('d:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analyses\spikeCode'));


for igroup=1:2
    if igroup==1
        rats={'LB01','LB03','LB05'};
    elseif igroup==2
        rats={'LB02','LB04','LB06'};
    end
    for irats=1:length(rats)
        parent=strcat('D:\ClarkP30_Recordings\',rats(irats));
        disp(['CYCLING THROUGH RAT:',char(rats(irats))])
        parent=char(parent);
        structdir=dir(parent);
        K=1;
        for I=1:length(structdir) % 1 TO # OF FILES IN DIR
            if structdir(I).isdir && structdir(I).name(1) ~= '.' % IF DIR NOT '.'
                if exist(([parent filesep structdir(I).name filesep 'TT']),'file'); % LOCATE TT FOLDER
                    CurentFolder=cd([parent filesep structdir(I).name filesep 'TT']); % CD TO TT FOLDER
                    CurrentMat=dir('*.mat'); % LOCATE .MAT FILES IN TT FOLDER
                    for J=1:length(CurrentMat) % 1 TO # OF .MAT FILES IN TT FOLDER
                        CurrentMatworking=CurrentMat(J).name;
                        if sum(ismember(CurrentMatworking,'spikeData'))==11 && sum(ismember(CurrentMatworking,'pathProperties'))~=16
                            warning('off','MATLAB:load:variableNotFound');
                            load(CurrentMatworking,'mean_vector_length','peak_Firing_Rate','BinsNbSpikes',...
                                'preferred_Direction','nSpikes','BinsAngle3','BinsAngle'); % OPEN .MAT WITH THESE VARS
                            disp(['READING: ',CurentFolder filesep CurrentMatworking])
                            if exist('mean_vector_length','var')
                                if mean_vector_length>=.5 && nSpikes>=50 % Arbitrary r length - will update with shuffled data once finished lb 2017
                                    
                                    rmpath('d:\Users\BClarkLab\Google Drive\MatlabDir\BClarkToolbox\Analyses\spikeCode');
                                    %using robust local regression 2nd degree polynomial,span based off 20% of data 
                                    smooth_BinsNbSpikes=smooth(BinsNbSpikes,.1,'rloess'); 
                                    addpath(genpath('d:\Users\BClarkLab\Google Drive\MatlabDir\BClarkToolbox\Analyses\spikeCode'))
                                    BinsNbSpikes(K,:)=
                                    BinsAngle3(K,:)=
                                    BinsAngle(K,:)=
                                    smooth_BinsNbSpikes(K,:)=                                    
                                    HD_cells(K,:)=[mean_vector_length,peak_Firing_Rate,preferred_Direction,nSpikes];
                                    K=K+1;
                                end
                            end
                            if ~exist('mean_vector_length','var') % IF THE .MAT FILE WAS INCORRECT
                                keep('I','J','K','structdir','parent','CurrentMat','CurrentFolder','rats','irats','igroup');
                            else
                                keep('I','J','K','HD_cells','structdir','parent','CurrentMat','CurentFolder','rats','irats','igroup'); % KEEP ALLSPIKEDATA & STRUCTDIR IF .MAT WAS CORRECT
                            end
                        end
                    end
                    keep('I','K','HD_cells','structdir','parent','CurrentMat','rats','irats','igroup'); % DON'T KEEP VAR J BECAUSE YOU JUST EXITED THE J LOOP
                end
            end
        end
        disp(['DONE WITH:',char(rats(irats))])
        keep('rats','irats','parent','igroup','HD_cells')
    end
    if igroup==1
        Group=('Control');
    elseif igroup==2
        Group=('Tg+');
    end
    
    % Plot Smoothed firing rate x HEAD DIRECTION
    figure (i), plot(BinsAngle3,smooth_BinsNbSpikes,'LineWidth',2,'color','k') %plot(BinsAngle3,BinsNbSpikes,'LineWidth',2,'color','k')
    axis tight;
    hold on
    xlim ([0 360]);
    box off
    fig3=figure (i+2);
    
    % Firing rate x HD polar plot for the nonsmoothed data above
    figure(i+1)
    polarplot = polar(BinsAngle([1:60 1]),BinsNbSpikes([1:60 1]),'b');
    set(polarplot, 'linewidth',3,'color','k');
    axis off
    title(['Polor Plot /',' MeanVecLength: ',num2str(mean_vector_length)]);
    set(0,'Showhiddenhandles','on')
    
    
    [filepath, filename] = fileparts(parent);
    saveas(fig,[filepath filesep Group,'_SmoothedRatePlot.tiff']);
    save([filepath filesep Group '_Data.mat']);
    keep('igroup')
    close all
end



