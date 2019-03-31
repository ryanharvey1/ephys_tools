% Saveallfigures
close all
cd D:\Projects\PAE_PlaceCell\ProcessedData
sessions=dir('*.mat');
sessions={sessions.name};

for i=1:length(sessions)
data=load(sessions{i}); 
p=postprocessFigures(data);

fig = get(groot,'CurrentFigure');

while ~isempty(ishandle(fig) & strcmp(get(fig, 'type'), 'figure'))

    id=strsplit(fig.Name,'  ');
    
    if exist(['D:\Projects\PAE_PlaceCell\Figures',filesep,id{1},filesep,id{2}],'dir')==0
        mkdir(['D:\Projects\PAE_PlaceCell\Figures',filesep,id{1},filesep,id{2}]);
    end
    
    set(fig, 'Position', get(0, 'Screensize'));
    if length(id)==3
        print(fig,'-dpng', '-r150',['D:\Projects\PAE_PlaceCell\Figures',filesep,id{1},filesep,id{2},filesep,erase(id{3},':'),'.png'])
    elseif length(id)==2
        print(fig,'-dpng', '-r150',['D:\Projects\PAE_PlaceCell\Figures',filesep,id{1},filesep,id{2},filesep,'raster.png'])
    end
    close
    fig = get(groot,'CurrentFigure');
end
end