% deeplabcuts_postprocess RH May 2019, edits by LB June 2019

figs=0;

%Initialize path to DLC output
path_to_files='d:\Users\BClarkLab\Desktop\Videos\lgOF_dlc';
files = dir([path_to_files,'\**\*.csv']);

%Initialize data table and variables
params=table;
params.subID{1}=[];
params.nose{1}=[];
params.head{1}=[];
params.back{1}=[];
params.butt{1}=[];
%Get SubID from subdirectory names
vidfile=dir(path_to_files);
vidfile(1:2,:)=[]; %get rid of . & .. entries
dia=202; %diameter of maze in cm

%% Compile environment Max/Min for transformation into cm and cue Coords for
%computing cue related measures in OF_postprocess

load('D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data\maxmin.mat');

for i=1:length(files)
    params.subID{i}=vidfile(i).name;
    params.dia{i}=dia;
    % load header
    fileID = fopen(fullfile(files(i).folder,files(i).name),'r');
    dataArray = textscan(fileID, '%q%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]',...
        3-2+1, 'Delimiter',',', 'TextType', 'string', 'HeaderLines',...
        2-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    header = [dataArray{1:end-1}];
    clearvars fileID dataArray ans;
    
    % load data
    fileID = fopen(fullfile(files(i).folder,files(i).name),'r');
    dataArray = textscan(fileID, '%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]',...
        'Delimiter', ',','TextType', 'string', 'EmptyValue', NaN,...
        'HeaderLines' ,4-1,'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    tsxy= [dataArray{1:end-1}];
    clearvars fileID dataArray ans;
    
    % load coded behaviors(obtained from Behavior_code.m)
    if exist(['G:\Maria\OF\Videos\lgOF_dlc\',vidfile(i).name,filesep,vidfile(i).name,'.mat'])
        behav=load(['G:\Maria\OF\Videos\lgOF_dlc\',vidfile(i).name,filesep,vidfile(i).name,'.mat']);
    else
        behav=nan;
    end
    
    % convert ts into seconds
    tsxy(:,1)=0:1/30:tsxy(end,1)/30;
    
    xunit=sin(0:pi/360:2*pi)*((maxmin.xmax{i}-maxmin.xmin{i})/2)*1.1+median([maxmin.xmax{i},maxmin.xmin{i}]);
    yunit=cos(0:pi/360:2*pi)*((maxmin.ymax{i}-maxmin.ymin{i})/2)*1.1+median([maxmin.ymax{i},maxmin.ymin{i}]);

    % Turn low liklihood point into NAN and smooth remaining coords
    idx=tsxy(:,contains(header(2,:),'likelihood'))<.95;
    xloc=find(contains(header(2,:),'x'));
    yloc=find(contains(header(2,:),'y'));
    tstart=behav.Trial_start(1,1);
    tend=tstart+(30*60); %30=30min trial duration, 60=seconds/min, 
    for l=1:size(idx,2)
        tsxy(idx(:,l),xloc(l):yloc(l))=NaN;
        
        tsxy(~inpolygon(tsxy(:,xloc(l)),tsxy(:,yloc(l)),xunit,yunit),xloc(l):yloc(l))=NaN;
        
        [tsxy(:,xloc(l)),tsxy(:,yloc(l))]=FixPos(tsxy(:,xloc(l)),tsxy(:,yloc(l)),tsxy(:,1));
    end
    
    %Keep only points that are within the trial start time plus 30min
    tsxy=tsxy(tsxy(:,1)>=tstart & tsxy(:,1)<=tend,:);
    
    %Reset timestamps to zero 
    tsxy(:,1)=(tsxy(:,1)-tsxy(1,1));
    
    % Save data
    params.ts{i}=tsxy(:,1);
    params.nose{i}(:,1) = tsxy(:,2);
    params.nose{i}(:,2) = tsxy(:,3);
    params.head{i}(:,1)=tsxy(:,5);
    params.head{i}(:,2)=tsxy(:,6);
    params.back{i}(:,1) = tsxy(:,8);
    params.back{i}(:,2) = tsxy(:,9);
    params.butt{i}(:,1)=tsxy(:,11);
    params.butt{i}(:,2)=tsxy(:,12);
    
    if figs==1
        figure;
        plot(params.nose{i}(:,1),params.nose{i}(:,2),'.k')
        hold on
        nPlot = plot(NaN,NaN,'ro','MarkerFaceColor','r');
        hPlot = plot(NaN,NaN,'bo','MarkerFaceColor','b');
        bPlot = plot(NaN,NaN,'go','MarkerFaceColor','g');
        buttPlot = plot(NaN,NaN,'co','MarkerFaceColor','c');
        
        
        xlim([min(params.nose{i}(:,1)) max(params.nose{i}(:,1))]);
        ylim([min(params.nose{i}(:,2)) max(params.nose{i}(:,2))]);
        % iterate through each point on line
        for k=1:length(params.nose{i}(:,1))
            set(nPlot,'XData',params.nose{i}(k,1),'YData',params.nose{i}(k,2));
            set(hPlot,'XData',params.head{i}(k,1),'YData',params.head{i}(k,2));
            set(bPlot,'XData',params.back{i}(k,1),'YData',params.back{i}(k,2));
            set(buttPlot,'XData',params.butt{i}(k,1),'YData',params.butt{i}(k,2));
            
            pause(1/240);
        end
        close
    end
end

load('D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data\cueCoords.mat');

params.Xmax=maxmin.xmax;
params.Xmin=maxmin.xmin;
params.Ymax=maxmin.ymax;
params.Ymin=maxmin.ymin;
params.cueCoords=cueCoords.coords;


