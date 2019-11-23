function visualize_cells(groupid,savelocation,varargin)
% visualize_cells: plots and saves all cell in groupid list
%
%   Input:  
%       groupid: cell array with same format as below (postprocessed .mat file,
%                   tetrode id, & cell number)
%
%     {'RH13_S20160808101427.mat','TT2.mat','2';
%     'RH13_S20160808101427.mat','TT2.mat','4';
%     'RH13_S20160808101427.mat','TT4.mat','5';
%     'RH13_S20160808101427.mat','TT5.mat','2';
%     'RH13_S20160808101427.mat','TT8.mat','3';
%     'RH13_S20160808103145.mat','TT4.mat','1'}
%       ...
%       savelocation: location with you want figure to be saved
%
%       Optional:
%                dpi: image quality in dots per inch (default 300)
%                format: image format (default png)
%   
%   Helper: add your processedData folder to path before running
%
% Ryan Harvey (2019)

p = inputParser;
p.addParameter('dpi','-r300');
p.addParameter('format','png');
p.parse(varargin{:});
dpi = p.Results.dpi;
format = p.Results.format;

if ~exist(savelocation,'file')
    mkdir(savelocation)
end

% find unique session: only loading sessions once significantly reduces run time 

sessions=unique(groupid(:,1));

for i=1:length(sessions)
    close all 
    idx=find(contains(groupid(:,1),sessions{i}));
    
    data=load(sessions{i});
    
    postprocessFigures.main(data,'cellid',{groupid(idx,2),str2double(groupid(idx,3))});
    
    FigList=findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig=1:length(FigList)
        FigHandle=FigList(iFig);
%         FigName=get(FigHandle, 'Name');
        
        set(FigHandle, 'Position', get(0, 'Screensize'));
        
        filename = FigHandle.Name;
        filename(isspace(filename)) = [];
        filename = strrep(filename, ':', '');
        filename = strrep(filename, '.ntt', '');

        print(FigHandle,['-d',format], dpi,...
            [savelocation,filesep,filename,'.',format])
    end
end
end