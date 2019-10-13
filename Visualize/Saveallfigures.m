% Saveallfigures
close all
% cd F:\Projects\PAE_PlaceCell\ProcessedData
cd ..
processedpath = pwd;
% processedpath=strsplit(data.session_path,filesep);
% processedpath(end-2:end)=[];
% cd(fullfile(strjoin(processedpath,filesep),'Figures'))
cd(fullfile(pwd,'ProcessedData'))
sessions=dir(fullfile(processedpath,'ProcessedData','*.mat'));
sessions={sessions.name};

for i=1:length(sessions)
    sess_id = extractBetween(sessions{i},'_','.mat');

    if exist(fullfile(processedpath,'Figures',extractBefore(sessions{i},'_'),sess_id{:}),'dir')
        continue
    end
    
    data=load(sessions{i});
    postprocessFigures.main(data,'colorcode','HD');
    
    fig = get(groot,'CurrentFigure');
    
    while ~isempty(ishandle(fig) & strcmp(get(fig, 'type'), 'figure'))
        
        id=strsplit(fig.Name,'  ');
        
        if exist([fullfile(processedpath,'Figures'),filesep,id{1},filesep,id{2}],'dir')==0
            mkdir([fullfile(processedpath,'Figures'),filesep,id{1},filesep,id{2}]);
        end
        
        set(fig, 'Position', get(0, 'Screensize'));
        if length(id)==3
            print(fig,'-dpng', '-r300',[fullfile(processedpath,'Figures'),filesep,id{1},filesep,id{2},filesep,erase(id{3},':'),'.png'])
        elseif length(id)==2
            print(fig,'-dpng', '-r300',[fullfile(processedpath,'Figures'),filesep,id{1},filesep,id{2},filesep,'raster.png'])
        end
        close
        fig = get(groot,'CurrentFigure');
    end
end