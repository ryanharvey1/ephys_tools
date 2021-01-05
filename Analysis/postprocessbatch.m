function postprocessbatch(dataset,varargin)
tic

p = inputParser;
p.addParameter('overwrite_lfp',0); % raw sample rate
p.addParameter('after_spike_sort_cleanup',0); % raw sample rate
p.parse(varargin{:});
overwrite_lfp = p.Results.overwrite_lfp;
after_spike_sort_cleanup = p.Results.after_spike_sort_cleanup;

if overwrite_lfp
    disp('WARNING. You are about to overwrite all neuroscope .xml and .lfp files')
    disp('Hit enter if this is okay')
    pause
end

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

close all
WaitMessage = parfor_wait(sum(nsessions),'Waitbar', true,'ReportInterval', 1);
% CYCLE THROUGH ALL DATA
parfor iratID=1:length(ratID)
    parent = strcat(dataset,filesep,ratID(iratID)); % CHANGE BASED ON PARENT FOLDER
    disp(['CYCLING THROUGH RAT:',char(ratID(iratID))])
    parent=char(parent);
    structdir=dir(parent);
    for I=1:length(structdir) % 1 TO # OF FILES IN DIR
        if structdir(I).isdir && structdir(I).name(1) ~= '.'  % IF NOT '.'
            basepath=[parent filesep structdir(I).name]; % SET PATH TO DATA THAT HAS BEEN THROUGH MCLUST
            if sum(ismember(sessions,['S',strjoin(regexp(structdir(I).name,'\d*','Match'),'')]))<1 &&...
                    exist([basepath,filesep,'Sorted'],'dir')==7
                
                % CALL after_spikesort_cleanup
                if after_spike_sort_cleanup
                    disp('running after_spikesort_cleanup')
                    delete(fullfile(basepath,'sorted\*.mat'))
                    after_spikesort_cleanup.main('path_name',basepath)
                end
                
                % CALL POSTPROCESS FUNCTION
                postprocess('path',basepath,'figures',0,'overwrite_lfp',overwrite_lfp); 
                
                % update wait 
                WaitMessage.Send;
                
            end
        end
    end
end
WaitMessage.Destroy;
end