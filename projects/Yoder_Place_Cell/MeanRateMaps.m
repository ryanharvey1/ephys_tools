% MeanRateMaps
clear ; close all ; clc
FigureLocation='/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/Figures_For_Paper';
% paths=importdata('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/GW- Field_StatsPlace Cells_Tilted_Mice_1_1.xlsx');
load('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/GW_Field_StatsPlaceCells_Tilted_Mice.mat')

path_group=data.textdata.Field_StatsPlaceCells0x2DTilted(:,1:2);
path_group(1,:)=[];
group=path_group(:,2);
control_paths=path_group(strcmp(group,'Control'),1);
tilted_paths=path_group(strcmp(group,'Tilted'),1);

% just first session
Sess1con=control_paths(1:5:length(control_paths));
Sess1tilt=tilted_paths(1:5:length(tilted_paths));

%% build 3d matrix of ratemaps
% for i=1:length(Sess1con)
%     % load in mat file
%     split=strsplit(Sess1con{i},'.');
%     matfiles=strcat(split(1),'.',split(2),'_Data.mat');
%     load(matfiles{1},'RateMap')
%     % check dims
%     [r,c]=size(RateMap);
%     if r>54; RateMap(r-54,:)=[]; elseif r<54; RateMap=[RateMap;NaN(54-r,c)]; end
%     if c>41; RateMap(:,c-41)=[]; elseif c<41; RateMap=[NaN(r,41-c),RateMap]; end
%     % add to 3d matrix
%     RateMapsCon(:,:,i)=RateMap;
% end
%
% for i=1:length(Sess1tilt)
%     % load in mat file
%     split=strsplit(Sess1tilt{i},'.');
%     matfiles=strcat(split(1),'.',split(2),'_Data.mat');
%     load(matfiles{1},'RateMap')
%     % check dims
%     [r,c]=size(RateMap);
%     if r>54;
%         RateMap(1:r-54,:)=[];
%     elseif r<54;
%         RateMap=[RateMap;NaN(54-r,c)];
%     end
%     [r,c]=size(RateMap);
%     if c>41;
%         RateMap(:,1:c-41)=[];
%     elseif c<41;
%         RateMap=[NaN(r,41-c),RateMap];
%     end
%     % add to 3d matrix
%     RateMapsTilt(:,:,i)=RateMap;
% end


%% LOCATE CENTER OF EACH MAP
for j=1:5
Sess1con=control_paths(j:5:length(control_paths));
boundFig=figure(1);
for i=1:length(Sess1con)
    % load in mat file
    split=strsplit(Sess1con{i},'.');
    matfiles=strcat(split(1),'.',split(2),'_Data.mat');
    if exist(matfiles{1},'file')~=2;continue;end
    load(matfiles{1},'field')
    % find field
    [row,col]=find(field);
    %find boundary
    k=boundary(row,col);
    bound(i).control=[col(k),row(k)];
    Ec(i,1,j)=Eccentricity(col(k),row(k));
%     boundFig;subplot(1,2,1);plot(col(k),row(k));hold on
    % find center point
    xmax=max(col);xmin=min(col);
    ymax=max(row);ymin=min(row);
    xc(i)=round(median(xmin:xmax));
    yc(i)=round(median(ymin:ymax));
%     subplot(1,2,1);s=scatter(xc(i),yc(i),'filled','r');hold on
%     s.LineWidth=30;
%     s.SizeData=1000;
end
end

for j=1:5
Sess1tilt=tilted_paths(j:5:length(tilted_paths));
for i=1:length(Sess1tilt)
    % load in mat file
    split=strsplit(Sess1tilt{i},'.');
    matfiles=strcat(split(1),'.',split(2),'_Data.mat');
    if exist(matfiles{1},'file')~=2;continue;end
    load(matfiles{1},'field')
    % find field
    [row,col]=find(field);
    %find boundary
    k=boundary(row,col);
    bound(i).tilted=[col(k),row(k)];
    Et(i,1,j)=Eccentricity(col(k),row(k));
%     boundFig;subplot(1,2,2);plot(col(k),row(k));hold on;
    % find center point
    xmax=max(col);xmin=min(col);
    ymax=max(row);ymin=min(row);
    xt(i)=round(median(xmin:xmax));
    yt(i)=round(median(ymin:ymax));
%     subplot(1,2,2);s=scatter(xt(i),yt(i),'filled','r');hold on
%     s.LineWidth=30;
%     s.SizeData=1000;
end
end
% figure;scatter(xc,yc,'*r');figure;scatter(xt,yt,'*r')
PlotsStat(Ec,Et,{'Eccentricity'})
% 
% x_max=max(xt)
% x_min=min(xt)
% y_max=max(yt)
% y_min=min(yt)
%             x=median([x_min:x_max]); y=median([y_min:y_max]); rad=(max(x_min-x_max,y_min-y_max))/2;
%             th = 0:pi/179.5:2*pi; % 0 to 2*pi(6.28318530717959) at 0.0175 increments to equal 360 points
%             xunit = rad * cos(th) + x;
%             yunit = rad * sin(th) + y;

%% 3D MATRIX OF FIELDS
for i=1:length(Sess1con)
    % load in mat file
    split=strsplit(Sess1con{i},'.');
    matfiles=strcat(split(1),'.',split(2),'_Data.mat');
    load(matfiles{1},'field')
    if ~exist('field','var');continue;end
    % check dims
    [r,c]=size(field);
    if r>28
        field(end-((r-28)-1):end,:)=[];
    elseif r<28
        field=[field;NaN(28-r,c)];
    end
    
    [r,c]=size(field);
    if c>26
        field(:,1:c-26)=[];
    elseif c<26
        field=[NaN(r,26-c),field];
    end
    % add to 3d matrix
    RateMapsCon(:,:,i)=field;
    clear field
end

for i=1:length(Sess1tilt)
    % load in mat file
    split=strsplit(Sess1tilt{i},'.');
    matfiles=strcat(split(1),'.',split(2),'_Data.mat');
    load(matfiles{1},'field')
    if ~exist('field','var');continue;end
    % check dims
    [r,c]=size(field);
    if r>28
        field(end-((r-28)-1):end,:)=[];
    elseif r<28
        field=[field;NaN(28-r,c)];
    end
    
    [r,c]=size(field);
    if c>26
        field(:,1:c-26)=[];
    elseif c<26
        field=[NaN(r,26-c),field];
    end
    % add to 3d matrix
    RateMapsTilt(:,:,i)=field;
    clear field
end
control=nansum(RateMapsCon,3);
% filtWidth = [3 3]; filtSigma = 1;
% imageFilter=fspecial('gaussian',filtWidth,filtSigma);
% control = nanconv(control,imageFilter, 'nanout');

tilted=nansum(RateMapsTilt,3);
% filtWidth = [3 3]; filtSigma = 1;
% imageFilter=fspecial('gaussian',filtWidth,filtSigma);
% tilted = nanconv(tilted,imageFilter, 'nanout');

%% LOCATE BORDERS
% [rows,columns,depth]=size(RateMapsTilt);
% for d=1:depth
%     for i=1:rows
%         for J=1:columns
%             if isnan(RateMapsTilt(i,J,d))==0
%                 index1(i)=i;
%             end
%         end
%     end
%     index1=index1(index1~=0); y_min(d)=index1(1); y_max(d)=numel(RateMapsTilt(:,1,d));
%
%     for i=1:columns
%         for J=1:rows
%             if isnan(RateMapsTilt(J,i,d))==0
%                 index(i)=i;
%             end
%         end
%     end
%     index=index(index~=0); x_min(d)=index(1); x_max(d)=numel(RateMapsTilt(1,:,d));
% end






meanRateMapcontrol=nanmean(RateMapsCon,3);
meanRateMaptilted=nanmean(RateMapsTilt,3);

meanRateMaptilted(meanRateMaptilted==0)=NaN;

alteredtiltedmap=meanRateMaptilted;

figure; h = pcolor(alteredtiltedmap);
colormap(jet);
axis square tight
hold on
set(h, 'EdgeColor', 'none');
box off
set(gca,'YDir','reverse');


% LOCATE BORDERS
[rows,columns]=size(meanRateMapcontrol);
for i=1:rows
    for J=1:columns
        if isnan(meanRateMapcontrol(i,J))==0
            index1(i)=i;
        end
    end
end
index1=index1(index1~=0); y_max=index1(1); y_min = numel(meanRateMapcontrol(:,1));

for i=1:columns
    for J=1:rows
        if isnan(meanRateMapcontrol(J,i))==0
            index(i)=i;
        end
    end
end
index=index(index~=0); x_max=index(1); x_min = numel(meanRateMapcontrol(1,:));

control=meanRateMapcontrol(y_max:y_min,x_max:x_min);
tilted=meanRateMaptilted(y_max:y_min,x_max:x_min);

% tilted(tilted>1)=tilted(tilted>1)+.5;

filtWidth = [3 3]; filtSigma = 1;
imageFilter=fspecial('gaussian',filtWidth,filtSigma);
control = nanconv(control,imageFilter, 'nanout');

filtWidth = [3 3]; filtSigma = 1;
imageFilter=fspecial('gaussian',filtWidth,filtSigma);
tilted = nanconv(tilted,imageFilter, 'nanout');

%%
control  = circRateMap(control);
tilted  = circRateMap(tilted);
%%


Average_Map=figure; Average_Map.Color=[1 1 1];
subplot(1,2,1);
h1 = pcolor(control);
colormap jet
% c=colorbar ;
% c.Visible='off';
axis square tight
hold on
set(h1, 'EdgeColor', 'none');
box off
set(gca,'YDir','reverse');
shading interp
ax=gca; ax.Visible='off';

subplot(1,2,2)
h2 = pcolor(tilted);
colormap jet
% c=colorbar ;
axis square tight
hold on
set(h2, 'EdgeColor', 'none');
box off
set(gca,'YDir','reverse');
shading interp
ax=gca; ax.Visible='off';

%%
print(Average_Map, '-dpng', '-r400',[FigureLocation,filesep,'Summed_Map.png'])

