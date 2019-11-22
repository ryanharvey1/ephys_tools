function postprocessbatch_rev(dataset)
tic

% check to see if the temp file exist...if not make it
cd (dataset)
cd ..

if exist([cd,filesep,'ProcessedData'],'dir')==0
    mkdir([cd,filesep,'ProcessedData'])
end

sessions=struct2table(dir( [cd,filesep,'ProcessedData',filesep,'*.mat']));
sessions=table2cell(sessions(:,1));
sessions=extractBetween(sessions,'_','.mat');

ratID=table2cell(struct2table(dir(dataset)));
paths=ratID(~contains([ratID(:,1)],'.'), 2);
ratID=ratID(~contains([ratID(:,1)],'.'), 1);
paths=strcat(paths,filesep,ratID);


for i=1:length(ratID)
    s=table2cell(struct2table(dir([dataset,filesep,ratID{i}])));
    nsessions(i)=sum(~contains([s(:,1)],'.'));
end
nsessions=sum(nsessions);ns=1;

close all
% h = waitbar(0,'Initializing waitbar...');

% flip rat list
ratID = flip(ratID);

% CYCLE THROUGH ALL DATA
for iratID=1:length(ratID)
    parent = strcat(dataset,filesep,ratID(iratID)); % CHANGE BASED ON PARENT FOLDER
    disp(['CYCLING THROUGH RAT:',char(ratID(iratID))])
    parent=char(parent);
    structdir=dir(parent);
    for I=1:length(structdir) % 1 TO # OF FILES IN DIR
        if structdir(I).isdir && structdir(I).name(1) ~= '.'  % IF NOT '.'
            path=[parent filesep structdir(I).name]; % SET PATH TO DATA THAT HAS BEEN THROUGH MCLUST
            if sum(ismember(sessions,['S',strjoin(regexp(structdir(I).name,'\d*','Match'),'')]))<1 && exist([path,filesep,'Sorted'],'dir')==7
                postprocess('path',path,'figures',0); % CALL POSTPROCESS FUNCTION
            end
%             waitbar(ns/nsessions,h,sprintf('%d%% done...',round((ns/nsessions)*100)))
            ns=ns+1;
        end
    end
    %     perc=(iratID/length(ratID))*100;
%     waitbar(perc/100,h,sprintf('%d%% done...',round(perc)))
end
% close(h)
toc
end