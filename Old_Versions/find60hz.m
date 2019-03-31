% find60hz
% This script cycles through all data and finds potential 60hz clusters
% Note... these are only potential clusters. Files this script locates might 
% or might not be 60hz. Further visual analysis must take place. 
%
% Ryan Harvey 2/1/2017
clear,clc
rats={'RH11','RH13','RH14','RH16'};
for irats=1:length(rats)
parent = strcat('D:\Place_Cell_Data\PAE_Rat\',rats(irats));
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
            CurentFolder=char(strcat(foldername(5)));
            CurrentMat=dir('*.mat'); % LOCATE .MAT FILES IN TT FOLDER
            for J=1:length(CurrentMat) % 1 TO # OF .MAT FILES IN TT FOLDER
                CurrentMatworking=CurrentMat(J).name; 
                try
                    load(CurrentMatworking,'InformationContent')
                    if exist('InformationContent','var')==1
                        [TTclust]= strsplit(CurrentMatworking,'_');
                        load(char(strcat(TTclust(1), '_ClusterSummary_', TTclust(2),'.mat')),'CI');
                        herts=nnz(CI.AutoCorr);
                        SpaceBetweenPeaks=mean(diff(find(CI.AutoCorr~=0)));                       
                        if herts>40 && herts<80 && SpaceBetweenPeaks>3.5 && SpaceBetweenPeaks<4.5 
                            CutupDir=strsplit(CI.CreatedInDirectory,filesep);
                            AllSpikeData(K,:)=strcat(strjoin(CutupDir(1:6),filesep),filesep, TTclust(1), '_ClusterSummary_', TTclust(2),'.mat'); % COMBINE INTO CELL ARRAY
                            K=K+1; % CYCLE THROUGH OUTPUT FILE SO YOU DON'T OVERWRITE
                        end
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
AllSpikeData = cell2table(AllSpikeData,'VariableNames',{'FolderWith60hz'}); % CREATE TABLE
[filepath, filename] = fileparts(parent);
if ismac==1
    writetable(AllSpikeData,[filepath filesep '_60hzPaths' filename '.csv'],'WriteRowNames',true);
else
    writetable(AllSpikeData,[filepath filesep '_60hzPaths' filename '.xlsx']);
end
disp(['DONE WITH:',char(rats(irats))])
keep('rats','irats','parent') 
end
disp('DONE')  