%% CompileAllMatFiles
% This script will cycle through a parent directory and pull out info from the .mat files
% Ryan Harvey 9/16
clear,clc
rats={'LB02'};
for irats=1:length(rats)
parent = strcat('G:\P30_Recordings\',rats(irats));
disp(['CYCLING THROUGH RAT:',char(rats(irats))])
addpath(genpath('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\Cell analysis')); 
parent=char(parent);
structdir=dir(parent);
K=1; 
for I=1:length(structdir) % 1 TO # OF FILES IN DIR
     if structdir(I).isdir && structdir(I).name(1) ~= '.' % IF DIR NOT '.'
        if exist(([parent filesep structdir(I).name filesep 'TT']),'file'); % LOCATE TT FOLDER
            CurentFolder=cd([parent filesep structdir(I).name filesep 'TT']); % CD TO TT FOLDER
            [foldername]= strsplit(CurentFolder,'\');
            CurentFolder=char(strcat(foldername(4)));
            CurrentMat=dir('*.mat'); % LOCATE .MAT FILES IN TT FOLDER
            for J=1:length(CurrentMat) % 1 TO # OF .MAT FILES IN TT FOLDER
                CurrentMatworking=CurrentMat(J).name; 
                try
                    load(CurrentMatworking,'InformationContent','sparsity','PeakRate','SpikeFile','spikeVec','DistanceFromTrackEnd','FieldWidth','mean_vector_length','preferred_Direction','Direct_infoContent'); % OPEN .MAT WITH THESE VARS
                    if exist('sparsity','var')==1
                        [TTclust]= strsplit(CurrentMatworking,'_');
                        load(char(strcat(TTclust(1), '_ClusterSummary_', TTclust(2),'.mat')),'CI');
                    end
                catch
                end % TRY CATCH END BASED ON DIFFERENT SPELLING AND OTHER MESSED UP SHIT
                try
                    load(CurrentMatworking,'InformationContent','Sparsity','PeakRate','SpikeFile','spikeVec','DistanceFromTrackEnd','FieldWidth','mean_vector_length','preferred_Direction','Direct_infoContent'); % FOR DIFFERENT SPELLINGS
                    if exist('Sparsity','var')==1
                        [TTclust]= strsplit(CurrentMatworking,'_');
                        load(char(strcat(TTclust(1), '_ClusterSummary_', TTclust(2),'.mat')),'CI');
                    end
                catch
                end
                try
                if exist('InformationContent','var')==1 % IF INFO CONTENT EXISTS IN A GIVEN .MAT FILE, IT'S PROBABLY THE RIGHT FILE
                   AllSpikeData(K,:)=[CurentFolder,CurrentMatworking,num2cell(InformationContent),num2cell(sparsity),num2cell(PeakRate),...
                   num2cell((length(SpikeFile)/length(spikeVec))*30),num2cell(max(DistanceFromTrackEnd)),num2cell(FieldWidth)...
                   num2cell(CI.PeakToTroughPts(1)),num2cell(CI.PeakToTroughPts(2)),num2cell(CI.PeakToTroughPts(3)),num2cell(CI.PeakToTroughPts(4)),...
                   num2cell(CI.PeakHalfWidthPts(1)),num2cell(CI.PeakHalfWidthPts(2)),num2cell(CI.PeakHalfWidthPts(3)),num2cell(CI.PeakHalfWidthPts(4)),...
                   num2cell(CI.nSpikes), num2cell(CI.CluSep.IsolationDist), num2cell(CI.CluSep.L_Ratio),num2cell(mean_vector_length),num2cell(preferred_Direction),num2cell(Direct_infoContent)]; % COMBINE INTO CELL ARRAY
                   K=K+1; % CYCLE THROUGH OUTPUT FILE SO YOU DON'T OVERWRITE
                end 
                catch
                end % TRY CATCH END BASED ON DIFFERENT SPELLING AND OTHER MESSED UP SHIT
                try
                if exist('InformationContent','var')==1 % IF INFO CONTENT EXISTS IN A GIVEN .MAT FILE, IT'S PROBABLY THE RIGHT FILE
                   AllSpikeData(K,:)=[CurentFolder,CurrentMatworking,num2cell(InformationContent),num2cell(Sparsity),num2cell(PeakRate),...
                   num2cell((length(SpikeFile)/length(spikeVec))*30),num2cell(max(DistanceFromTrackEnd)),num2cell(FieldWidth),...
                   num2cell(CI.PeakToTroughPts(1)),num2cell(CI.PeakToTroughPts(2)),num2cell(CI.PeakToTroughPts(3)),num2cell(CI.PeakToTroughPts(4)),...
                   num2cell(CI.PeakHalfWidthPts(1)),num2cell(CI.PeakHalfWidthPts(2)),num2cell(CI.PeakHalfWidthPts(3)),num2cell(CI.PeakHalfWidthPts(4)),...
                   num2cell(CI.nSpikes), num2cell(CI.CluSep.IsolationDist), num2cell(CI.CluSep.L_Ratio),num2cell(mean_vector_length),num2cell(preferred_Direction),num2cell(Direct_infoContent)]; % COMBINE INTO CELL ARRAY
                   K=K+1; % CYCLE THROUGH OUTPUT FILE SO YOU DON'T OVERWRITE
                end 
                catch
                end
                if ~exist('AllSpikeData','var') % IF THE .MAT FILE WAS INCORRECT
                    keep('I','J','K','structdir','parent','CurrentMat','CurrentFolder','rats','irats'); 
                else
                    keep('I','J','K','AllSpikeData','structdir','parent','CurrentMat','CurentFolder','rats','irats'); % KEEP ALLSPIKEDATA & STRUCTDIR IF .MAT WAS CORRECT
                end
            end
            keep('I','K','AllSpikeData','structdir','parent','CurrentMat','rats','irats'); % DON'T KEEP VAR J BECAUSE YOU JUST EXITED THE J LOOP
        end
     end
end
AllSpikeData = cell2table(AllSpikeData,'VariableNames',{'CurentFolder','SpikeFile','InformationContent','Sparsity','PeakRate','OverallFiringRate','DistanceFromTrackEnd','FieldWidth',...
    'PeaktoTrough1','PeaktoTrough2','PeaktoTrough3','PeaktoTrough4','PeakHalfWidthPts1','PeakHalfWidthPts2','PeakHalfWidthPts3','PeakHalfWidthPts4',...
    'nSpikes','IsolationDist','L_Ratio','mean_vector_length','preferred_Direction','Direct_infoContent'}); % CREATE TABLE
[filepath, filename] = fileparts(parent);
if ismac==1
    writetable(AllSpikeData,[filepath filesep '_AllSpikeData' filename '.csv'],'WriteRowNames',true);
else
    writetable(AllSpikeData,[filepath filesep '_AllSpikeData' filename '.xlsx']);
end
disp(['DONE WITH:',char(rats(irats))])
keep('rats','irats','parent') 
end

disp('done')    
   

