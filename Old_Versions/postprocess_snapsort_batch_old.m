% postprocess_snapsort_batch({'RH11','RH16','LS17','LS19','LE2813','RH13','RH14','LS21','LS23','LE2821','LE2823'})
%

%        Runs postprocess_snapsort script and saves the results to temp
%        files. It then loads and saves all the temp files back to a single data
%        stucutre 


function postprocess_snapsort_batch(rats)
addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis'));
addpath(genpath('d:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis'));

datafolder='D:\Place_Cell_Data\RawPAE_PlaceCell\';

% check to see if the temp file exist...if not make it
if exist('D:\Place_Cell_Data\RawPAE_PlaceCell\temp','dir')==0
    mkdir('D:\Place_Cell_Data\RawPAE_PlaceCell\temp')
end

sessions=struct2table(dir( 'D:\Place_Cell_Data\RawPAE_PlaceCell\temp\*.mat'));
sessions=table2cell(sessions(:,1));
sessions=extractBetween(sessions,'_','.mat');


for i=1:length(rats)
    nsessions(i)=length(dir([datafolder,rats{i}]));
end
nsessions=sum(nsessions);ns=1;

close all
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
            cd 'F:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox'
            if sum(ismember(sessions,['S',strjoin(regexp(structdir(I).name,'\d*','Match'),'')]))<1 && exist([path,filesep,'SNAPSorterResults'],'dir')==7
                track_length=TrackLength(path); % SET TRACK LENGTH
                postprocess_snapsort(path,track_length,'yes',0); % CALL POSTPROCESS FUNCTION
            end
            waitbar(ns/nsessions,h,sprintf('%d%% done...',round((ns/nsessions)*100)))
            ns=ns+1;
        end
    end
    %     perc=(irats/length(rats))*100;
%     waitbar(perc/100,h,sprintf('%d%% done...',round(perc)))
end
close(h)

disp('DONE WITH POST PROCESSING...')
disp('............................')
disp('............................')
disp('............................')
disp('RECOMPILING RESULTS INTO SINGLE DATA FILE')

sessions=struct2table(dir( 'D:\Place_Cell_Data\RawPAE_PlaceCell\temp\*.mat'));
sessions=table2cell(sessions(:,1));

sessionID=extractBetween(sessions,'_','.mat');
ratID=extractBefore(sessions,'_');

cd D:\Place_Cell_Data\RawPAE_PlaceCell\temp
for i=1:length(sessions)
    disp(sessions{i})
    load(sessions{i})
    datatemp.(ratID{i}).(sessionID{i})=data.(ratID{i}).(sessionID{i});
end
data=datatemp;
save('D:\Place_Cell_Data\RawPAE_PlaceCell\data.mat','data','-v7.3')
end



% OLD VERSION \/ \/ \/ \/

% 
% function postprocess_snapsort_batch(rats)
% addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis'));
% addpath(genpath('d:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis'));
% 
% datafolder='D:\Place_Cell_Data\RawPAE_PlaceCell\';
% linear_track = 'yes';
% load('D:\Place_Cell_Data\RawPAE_PlaceCell\data.mat','data')
% ratID=fieldnames(data);
% sessions=[];
% for i=1:length(ratID)
%     sessions=[sessions;fieldnames(data.(ratID{i}))];
% end
% 
% h = waitbar(0,'Initializing waitbar...');
% 
% % CYCLE THROUGH ALL DATA
% for irats=1:length(rats)
%     parent = strcat(datafolder,rats(irats)); % CHANGE BASED ON PARENT FOLDER
%     disp(['CYCLING THROUGH RAT:',char(rats(irats))])
%     parent=char(parent);
%     structdir=dir(parent);
%     for I=1:length(structdir) % 1 TO # OF FILES IN DIR
%         if structdir(I).isdir && structdir(I).name(1) ~= '.'  % IF NOT '.'
%             path=[parent filesep structdir(I).name]; % SET PATH TO DATA THAT HAS BEEN THROUGH MCLUST
%             cd 'F:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox'
%             if sum(ismember(sessions,['S',strjoin(regexp(structdir(I).name,'\d*','Match'),'')]))<1 && exist([path,filesep,'SNAPSorterResults'],'dir')==7
%                 track_length=TrackLength(path); % SET TRACK LENGTH
%                 postprocess_snapsort(path,track_length,linear_track,0); % CALL POSTPROCESS FUNCTION
%             end
%         end
%     end
%     perc=(irats/length(rats))*100;
%     waitbar(perc/100,h,sprintf('%d%% done...',round(perc)))
% end
% close(h)
% end
% 
% 


