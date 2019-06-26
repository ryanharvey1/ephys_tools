
%%
%%%%%%%%%%% Plot Paths with home base outline and cue location
% load('params_V8');
load('params_V16');
param_idx=params.subID;

vars=fieldnames(params);

%% Results Fig. 1 (Path Length, Search Area, Time in outer wall, running speed)
[row,~]=find(ismember(vars,{'pathL','searchArea','runSpeed'}));

idx=1;
fig=figure; fig.Color=[1 1 1];

for i=1:size(row,1)
tg1=vertcat(params.(vars{row(i)}){contains(param_idx,'Tg') & contains(param_idx,'D1')});
wt1=vertcat(params.(vars{row(i)}){contains(param_idx,'WT') & contains(param_idx,'D1')});
tg2=vertcat(params.(vars{row(i)}){contains(param_idx,'Tg') & contains(param_idx,'D2')});
wt2=vertcat(params.(vars{row(i)}){contains(param_idx,'WT') & contains(param_idx,'D2')});

subplot(size(row,1),2,idx)
plotspread_wrapper(wt1,tg1,{'WT','Tg'})
ylim([0 max([wt2;tg2;wt1;tg1])])
ylabel(vars{row(i)})
title('Proximal Cue Present')

idx=idx+1;
subplot(size(row,1),2,idx)
plotspread_wrapper(wt2,tg2,{'WT','Tg'})
ylim([0 max([wt2;tg2;wt1;tg1])])
title('Proximal Cue Absent')
idx=idx+1;
end

%% Figure 2

tg1=median(vertcat(params.dwellQuad{contains(param_idx,'Tg') & contains(param_idx,'D1')}),1);
wt1=median(vertcat(params.dwellQuad{contains(param_idx,'WT') & contains(param_idx,'D1')}),1);
tg2=median(vertcat(params.dwellQuad{contains(param_idx,'Tg') & contains(param_idx,'D2')}),1);
wt2=median(vertcat(params.dwellQuad{contains(param_idx,'WT') & contains(param_idx,'D2')}),1);

fig=figure; fig.Color=[1 1 1];
subplot(2,2,1) 
imagesc(tg1); colormap(viridis(255))
title('TgF344-AD Cue Present')
xlabel('Maze Quadrants (22.5deg/bin)')
subplot(2,2,3) 
imagesc(wt1); colormap(viridis(255))
title('Wild Type Cue Present')
xlabel('Maze Quadrants (22.5deg/bin)')
subplot(2,2,2) 
imagesc(tg2); colormap(viridis(255))
title('TgF344-AD Cue Absent')
xlabel('Maze Quadrants (22.5deg/bin)')
subplot(2,2,4) 
imagesc(wt2); colormap(viridis(255))
title('Wild Type Cue Absent')
xlabel('Maze Quadrants (22.5deg/bin)')

%% Surface plots for occupancy 
tg1=mean(cat(3,params.rateMap{contains(param_idx,'Tg') & contains(param_idx,'D1')}),3);
wt1=mean(cat(3,params.rateMap{contains(param_idx,'WT') & contains(param_idx,'D1')}),3);
tg2=mean(cat(3,params.rateMap{contains(param_idx,'Tg') & contains(param_idx,'D2')}),3);
wt2=mean(cat(3,params.rateMap{contains(param_idx,'WT') & contains(param_idx,'D2')}),3);

%Use meshgrid to serve as basis for logical mask.
imageSizeX = size(tg1,1);
imageSizeY = size(tg1,2);
[columnsInImage,rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);

clear imageSizeX imageSizeY

% Next create the circle in the image.
centerX = median(1:size(tg1,1)); centerY = median(1:size(tg1,2)); radius = median(1:size(tg1,2));
circlePixels = (rowsInImage - centerY).^2 ...
    + (columnsInImage - centerX).^2 <= radius.^2;

clear rowsInImage columnsInImage

tg1(~circlePixels)=NaN; %indicate area outside the maze by labelling with NaN
wt1(~circlePixels)=NaN; %indicate area outside the maze by labelling with NaN
tg2(~circlePixels)=NaN; %indicate area outside the maze by labelling with NaN
wt2(~circlePixels)=NaN; %indicate area outside the maze by labelling with NaN

imAlpha1=ones(size(tg1));
imAlpha1(isnan(tg1))=0;
imAlpha2=ones(size(wt1));
imAlpha2(isnan(wt1))=0;
imAlpha3=ones(size(tg2));
imAlpha3(isnan(tg2))=0;
imAlpha4=ones(size(wt2));
imAlpha4(isnan(wt2))=0;

fig=figure; fig.Color=[1 1 1];
subaxis(2,2,1) 
surf(tg1,'AlphaData',imAlpha1); colormap(viridis(255));axis xy;hold on; box off; axis image;shading flat;
title('TgF344-AD Cue Present')
subaxis(2,2,3) 
surf(wt1,'AlphaData',imAlpha2); colormap(viridis(255));axis xy;hold on; box off; axis image;shading flat;
title('Wild Type Cue Present')
subaxis(2,2,2) 
surf(tg2,'AlphaData',imAlpha3); colormap(viridis(255));axis xy;hold on; box off; axis image;shading flat;
title('TgF344-AD Cue Absent')
subaxis(2,2,4) 
surf(wt2,'AlphaData',imAlpha4);colormap(viridis(255));axis xy;hold on; box off; axis image;shading flat;
title('Wild Type Cue Absent')

% 
% 
% fig=figure; fig.Color=[1 1 1];
% 
% subplot(2,2,1) 
% tg1(any(isnan(tg1)'),:)=[];
% scatter_COL=heatscatter(tg1(:,1),tg1(:,2));
% scatter(tg1(:,1),tg1(:,2), 150, smoothdata(scatter_COL,'gaussian',200), '.');
% colormap(viridis(10));axis xy; axis off; hold on; box off; axis image;
% title('TgF344-AD Cue Present')
% 
% subplot(2,2,3) 
% wt1(any(isnan(wt1)'),:)=[];
% scatter_COL=heatscatter(wt1(:,1),wt1(:,2));
% scatter(wt1(:,1),wt1(:,2), 150, smoothdata(scatter_COL,'gaussian',200), '.');
% colormap(viridis(10));axis xy; axis off; hold on; box off; axis image;
% title('Wild Type Cue Present')
% 
% subplot(2,2,2) 
% tg2(any(isnan(tg2)'),:)=[];
% scatter_COL=heatscatter(tg2(:,1),tg2(:,2));
% scatter(tg2(:,1),tg2(:,2), 150, smoothdata(scatter_COL,'gaussian',200), '.');
% colormap(viridis(10));axis xy; axis off; hold on; box off; axis image;
% title('TgF344-AD Cue Absent')
% 
% subplot(2,2,4) 
% wt2(any(isnan(wt2)'),:)=[];
% scatter_COL=heatscatter(wt2(:,1),wt2(:,2));
% scatter(wt2(:,1),wt2(:,2), 150, smoothdata(scatter_COL,'gaussian',200), '.');
% colormap(viridis(10));axis xy; axis off; hold on; box off; axis image;
% title('Wild Type Cue Absent')

%% Cross corr of occupancy maps
row=1;
for i=1:2:size(params,1)
    occCorr(row,1)=corr2(params.rateMap{i},params.rateMap{i+1});
    row=row+1;
end

figure; 
plotspread_wrapper(occCorr(13:end,1),occCorr(1:12,1),{'WT','Tg'})
title('Cross Corr between occupancy maps')

%% Proximity of primary hb day 1 versus primary hb day 2
row=1;
for i=1:2:size(params,1)
    primaryHBdist(row,1)=sqrt((params.HBcenter{i,1}{1,occIdx(i,1)}(1,1)-params.HBcenter{i+1,1}{1,occIdx(i+1,1)}(1,1))^2+...
        (params.HBcenter{i,1}{1,occIdx(i,1)}(1,2)-params.HBcenter{i+1,1}{1,occIdx(i+1,1)}(1,2))^2);
    row=row+1;
end

tgM=[1:5,14];
tgF=6:11;
wtM=[13:18];
wtF=[19:24];

figure; 
plotspread_wrapper(primaryHBdist(13:end,1),primaryHBdist(1:12,1),{'WT','Tg'})
title('Distance between Day 1 & Day 2 Primary Home Base Centroids')

figure;
scatter(occCorr(wtM,1),primaryHBdist(wtM,1),'k','filled')
hold on; 
scatter(occCorr(wtF,1),primaryHBdist(wtF,1),'MarkerEdgeColor',[.5 .5 .5],...
              'MarkerFaceColor',[.5 .5 .5])
hold on;
scatter(occCorr(tgM,1),primaryHBdist(tgM,1),'r','filled')
hold on; 
scatter(occCorr(tgF,1),primaryHBdist(tgF,1),'b','filled')
xlabel('Occupancy Map Correlation')
ylabel('Distace between Primary HB centers (D1vsD2)')
