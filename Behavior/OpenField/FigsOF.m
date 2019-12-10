
%%
%%%%%%%%%%% Plot Paths with home base outline and cue location
% load('params_V8');
load('D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data\params_V18');
param_idx=params.subID;

vars=fieldnames(params);

%%  Results Fig. 1 (Path Length, Search Area, Time in outer wall, running speed)
[row,~]=find(ismember(vars,{'pathL','searchArea','runSpeed','NumStops'}));
varlabel={'Path Length (cm)','Running Speed (cm/s)','Number of Stops','Search Area'};
idx=1;
fig=figure; fig.Color=[1 1 1];

for i=1:size(row,1)
tg1=vertcat(params.(vars{row(i)}){contains(param_idx,'Tg') & contains(param_idx,'D1')});
wt1=vertcat(params.(vars{row(i)}){contains(param_idx,'WT') & contains(param_idx,'D1')});
tg2=vertcat(params.(vars{row(i)}){contains(param_idx,'Tg') & contains(param_idx,'D2')});
wt2=vertcat(params.(vars{row(i)}){contains(param_idx,'WT') & contains(param_idx,'D2')});

subplot(size(row,1),2,idx)
plotspread_wrapper(wt1,tg1,{'F344','TgF344-AD'})
ylim([0 max([wt2;tg2;wt1;tg1])])
ylabel(varlabel{i})
title('Test')
idx=idx+1;
subplot(size(row,1),2,idx)
plotspread_wrapper(wt2,tg2,{'F344','TgF344-AD'})
ylim([0 max([wt2;tg2;wt1;tg1])])
title('Probe')
idx=idx+1;
end

export_fig('D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Figures\PNGs\whole_trial_meas.png','-m4') 

%%
for i=1:size(row,1)
tg1=vertcat(params.(vars{row(i)}){contains(param_idx,'Tg') & contains(param_idx,'D1')});
wt1=vertcat(params.(vars{row(i)}){contains(param_idx,'WT') & contains(param_idx,'D1')});
tg2=vertcat(params.(vars{row(i)}){contains(param_idx,'Tg') & contains(param_idx,'D2')});
wt2=vertcat(params.(vars{row(i)}){contains(param_idx,'WT') & contains(param_idx,'D2')});
figure; 
stat_plot(tg1,wt1,{'Tg','WT'},vars{row(i)},'plots',2)
stat_plot(tg2,wt2,{'Tg','WT'},vars{row(i)},'plots',2)

end

%%
tg1=horzcat(params.time2HB{contains(param_idx,'Tg') & contains(param_idx,'D1')});
wt1=horzcat(params.time2HB{contains(param_idx,'WT') & contains(param_idx,'D1')});
tg2=horzcat(params.time2HB{contains(param_idx,'Tg') & contains(param_idx,'D2')});
wt2=horzcat(params.time2HB{contains(param_idx,'WT') & contains(param_idx,'D2')});

fig = figure;
fig.Color = [1 1 1];
subplot(1,2,1)
[h1,h2]=CoolHistogram(cell2mat(wt1)',cell2mat(tg1)',10,{'Time to Home Base (s)'})
title('Test')
legend({'F344','TgF344-AD'})
ylim([0 1])
legend boxoff
subplot(1,2,2)
[h1,h2]=CoolHistogram(cell2mat(wt2)',cell2mat(tg2)',10,{'Time to Home Base (s)'})
title('Probe')
legend({'F344','TgF344-AD'})
legend boxoff
ylim([0 1])

%% Figure 2 - Median dwell time per qudrant
tg1=median(vertcat(params.dwellQuad{contains(param_idx,'Tg') & contains(param_idx,'D1')}),1);
wt1=median(vertcat(params.dwellQuad{contains(param_idx,'WT') & contains(param_idx,'D1')}),1);
tg2=median(vertcat(params.dwellQuad{contains(param_idx,'Tg') & contains(param_idx,'D2')}),1);
wt2=median(vertcat(params.dwellQuad{contains(param_idx,'WT') & contains(param_idx,'D2')}),1);

%Histograms for quadrant dwell time
f=figure; f.Color=[1 1 1];
subplot(1,2,1)
bar(sum(wt1,1)/sum(wt1(:)),'FaceColor','k','FaceAlpha',.5);
title('Cue Present')
ylim([0 1])
xlabel('Quadrant')
ylabel('Normalized Dwell Time')
 hold on 
 bar(sum(tg1,1)/sum(tg1(:)),'FaceColor','r','FaceAlpha',.5);
 set(gca,'FontSize',14,'FontWeight','bold','FontName','Calibri')
box off 

subplot(1,2,2)
bar(sum(wt2,1)/sum(wt2(:)),'FaceColor','k','FaceAlpha',.5);
title('Cue Absent')
ylim([0 1])
xlabel('Quadrant')
ylabel('Normalized Dwell Time')
 hold on 
 bar(sum(tg2,1)/sum(tg2(:)),'FaceColor','r','FaceAlpha',.5);
  set(gca,'FontSize',14,'FontWeight','bold','FontName','Calibri')
box off 


%% Bar graphs for Ben 

% Figure - Bar chart for quadrants - Dweel time per quadrant (sec)
tg1=mean(vertcat(params.dwellQuad{contains(param_idx,'Tg') & contains(param_idx,'D1')}),1);
wt1=mean(vertcat(params.dwellQuad{contains(param_idx,'WT') & contains(param_idx,'D1')}),1);
tg2=mean(vertcat(params.dwellQuad{contains(param_idx,'Tg') & contains(param_idx,'D2')}),1);
wt2=mean(vertcat(params.dwellQuad{contains(param_idx,'WT') & contains(param_idx,'D2')}),1);

f=figure; f.Color=[1 1 1];
subplot(5,2,1)
bar(wt1,'FaceColor','k','FaceAlpha',.5);
title('Cue Present')
xlabel('Quadrant')
ylabel('Dwell(s)')
 hold on 
 bar(tg1,'FaceColor','r','FaceAlpha',.5);
 set(gca,'FontSize',11,'FontWeight','bold','FontName','Calibri')
box off 

subplot(5,2,2)
bar(wt2,'FaceColor','k','FaceAlpha',.5);
title('Cue Absent')
xlabel('Quadrant')
ylabel('Dwell (s)')
 hold on 
 bar(tg2,'FaceColor','r','FaceAlpha',.5);
  set(gca,'FontSize',11,'FontWeight','bold','FontName','Calibri')
box off 

% Figure - Bar chart for quadrants - average angular velocity (theta/sec)
tg1=nanmean(abs(vertcat(params.angVelQuad{contains(param_idx,'Tg') & contains(param_idx,'D1')})),1);
wt1=nanmean(abs(vertcat(params.angVelQuad{contains(param_idx,'WT') & contains(param_idx,'D1')})),1);
tg2=nanmean(abs(vertcat(params.angVelQuad{contains(param_idx,'Tg') & contains(param_idx,'D2')})),1);
wt2=nanmean(abs(vertcat(params.angVelQuad{contains(param_idx,'WT') & contains(param_idx,'D2')})),1);

subplot(5,2,3)
bar(wt1,'FaceColor','k','FaceAlpha',.5);
title('Cue Present')
xlabel('Quadrant')
ylabel('Mean Ang Vel (deg/s)')
 hold on 
 bar(tg1,'FaceColor','r','FaceAlpha',.5);
 set(gca,'FontSize',11,'FontWeight','bold','FontName','Calibri')
box off 

subplot(5,2,4)
bar(wt2,'FaceColor','k','FaceAlpha',.5);
title('Cue Absent')
xlabel('Quadrant')
ylabel('Mean Ang Vel (deg/s)')
 hold on 
 bar(tg2,'FaceColor','r','FaceAlpha',.5);
  set(gca,'FontSize',11,'FontWeight','bold','FontName','Calibri')
box off 

% Figure - Bar chart for quadrants - Average running speed (cm/sec)
tg1=nanmean(vertcat(params.pathIVQuad{contains(param_idx,'Tg') & contains(param_idx,'D1')}),1);
wt1=nanmean(vertcat(params.pathIVQuad{contains(param_idx,'WT') & contains(param_idx,'D1')}),1);
tg2=nanmean(vertcat(params.pathIVQuad{contains(param_idx,'Tg') & contains(param_idx,'D2')}),1);
wt2=nanmean(vertcat(params.pathIVQuad{contains(param_idx,'WT') & contains(param_idx,'D2')}),1);


subplot(5,2,5)
bar(wt1,'FaceColor','k','FaceAlpha',.5);
title('Cue Present')
xlabel('Quadrant')
ylabel('Average lin vel (cm/s)')
 hold on 
 bar(tg1,'FaceColor','r','FaceAlpha',.5);
 set(gca,'FontSize',11,'FontWeight','bold','FontName','Calibri')
box off 

subplot(5,2,6)
bar(wt2,'FaceColor','k','FaceAlpha',.5);
title('Cue Absent')
xlabel('Quadrant')
ylabel('Average lin vel (cm/s)')
 hold on 
 bar(tg2,'FaceColor','r','FaceAlpha',.5);
  set(gca,'FontSize',11,'FontWeight','bold','FontName','Calibri')
box off 

% Figure - Bar chart for quadrants - number of stops per Quadrant
tg1=nanmean(vertcat(params.numstopQuad{contains(param_idx,'Tg') & contains(param_idx,'D1')}),1);
wt1=nanmean(vertcat(params.numstopQuad{contains(param_idx,'WT') & contains(param_idx,'D1')}),1);
tg2=nanmean(vertcat(params.numstopQuad{contains(param_idx,'Tg') & contains(param_idx,'D2')}),1);
wt2=nanmean(vertcat(params.numstopQuad{contains(param_idx,'WT') & contains(param_idx,'D2')}),1);


subplot(5,2,7)
bar(wt1,'FaceColor','k','FaceAlpha',.5);
title('Cue Present')
xlabel('Quadrant')
ylabel('Average Number of Stops')
 hold on 
 bar(tg1,'FaceColor','r','FaceAlpha',.5);
 set(gca,'FontSize',11,'FontWeight','bold','FontName','Calibri')
box off 

subplot(5,2,8)
bar(wt2,'FaceColor','k','FaceAlpha',.5);
title('Cue Absent')
xlabel('Quadrant')
ylabel('Average Number of Stops')
 hold on 
 bar(tg2,'FaceColor','r','FaceAlpha',.5);
  set(gca,'FontSize',11,'FontWeight','bold','FontName','Calibri')
box off 

% Figure - Bar chart for quadrants - Path Length (cm)
tg1=nanmean(vertcat(params.pathLQuad{contains(param_idx,'Tg') & contains(param_idx,'D1')}),1);
wt1=nanmean(vertcat(params.pathLQuad{contains(param_idx,'WT') & contains(param_idx,'D1')}),1);
tg2=nanmean(vertcat(params.pathLQuad{contains(param_idx,'Tg') & contains(param_idx,'D2')}),1);
wt2=nanmean(vertcat(params.pathLQuad{contains(param_idx,'WT') & contains(param_idx,'D2')}),1);


%Histograms for quadrant dwell time
subplot(1,2,1)
bar(wt1,'FaceColor','k','FaceAlpha',.5);
title('Cue Present')
xlabel('Quadrant')
ylabel('Path Length (cm)')
 hold on 
 bar(tg1,'FaceColor','r','FaceAlpha',.5);
 set(gca,'FontSize',11,'FontWeight','bold','FontName','Calibri')
box off 

subplot(1,2,2)
bar(wt2,'FaceColor','k','FaceAlpha',.5);
title('Cue Absent')
xlabel('Quadrant')
ylabel('Path Length (cm)')
 hold on 
 bar(tg2,'FaceColor','r','FaceAlpha',.5);
  set(gca,'FontSize',11,'FontWeight','bold','FontName','Calibri')
box off 




%% Colormap of quadrant dwell time - linearized

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
imagesc(tg1,'AlphaData',imAlpha1); colormap(viridis(255));axis xy;hold on; box off; axis image;shading flat;
title('TgF344-AD Test')
subaxis(2,2,3) 
imagesc(wt1,'AlphaData',imAlpha2); colormap(viridis(255));axis xy;hold on; box off; axis image;shading flat;
title('Wild Type Test')
subaxis(2,2,2) 
imagesc(tg2,'AlphaData',imAlpha3); colormap(viridis(255));axis xy;hold on; box off; axis image;shading flat;
title('TgF344-AD Probe')
subaxis(2,2,4) 
imagesc(wt2,'AlphaData',imAlpha4);colormap(viridis(255));axis xy;hold on; box off; axis image;shading flat;
title('Wild Type Probe')
axis off

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

figure; 
plotspread_wrapper(primaryHBdist(13:end,1),primaryHBdist(1:12,1),{'WT','Tg'})
title('Distance between Day 1 & Day 2 Primary Home Base Centroids')

figure; 
histogram(primaryHBdist(13:end,1),10,'FaceColor','k','FaceAlpha',.5); 
hold on; histogram(primaryHBdist(1:12,1),10,'FaceColor','r','FaceAlpha',.5)
legend({'F344','TgF344-AD'})


%% Stops by stop duration (size) and time in trial (color)
for i=1:size(params.tsStop,1) 
    for ii=1:size(params.tsStop{i},2)
        params.tsStopIdx{i}(1,ii)=params.tsStop{i}{1,ii}(1,1);
    end
end


tg1=vertcat(params.stopCenter{contains(param_idx,'Tg') & contains(param_idx,'D1')});
tg1_c=horzcat(params.timeStopped{contains(param_idx,'Tg') & contains(param_idx,'D1')});
tg1_ts=horzcat(params.tsStopIdx{contains(param_idx,'Tg') & contains(param_idx,'D1')});

wt1=vertcat(params.stopCenter{contains(param_idx,'WT') & contains(param_idx,'D1')});
wt1_c=horzcat(params.timeStopped{contains(param_idx,'WT') & contains(param_idx,'D1')});
wt1_ts=horzcat(params.tsStopIdx{contains(param_idx,'WT') & contains(param_idx,'D1')});

tg2=vertcat(params.stopCenter{contains(param_idx,'Tg') & contains(param_idx,'D2')});
tg2_c=horzcat(params.timeStopped{contains(param_idx,'Tg') & contains(param_idx,'D2')});
tg2_ts=horzcat(params.tsStopIdx{contains(param_idx,'Tg') & contains(param_idx,'D2')});

wt2=vertcat(params.stopCenter{contains(param_idx,'WT') & contains(param_idx,'D2')});
wt2_c=horzcat(params.timeStopped{contains(param_idx,'WT') & contains(param_idx,'D2')});
wt2_ts=horzcat(params.tsStopIdx{contains(param_idx,'WT') & contains(param_idx,'D2')});

fig=figure; 
fig.Color=[1 1 1];

subaxis(2,2,1)
plot(sin(0:2*pi/1000:2*pi)*101,cos(0:2*pi/1000:2*pi)*101,'k'); hold on;
scatter(tg1(:,1),tg1(:,2),cell2mat(tg1_c)',tg1_ts','Filled')
title('TgF344-AD Stops - Cue Present')
axis image 
axis off
colorbar
colormap(viridis(255))

subaxis(2,2,2)
plot(sin(0:2*pi/1000:2*pi)*101,cos(0:2*pi/1000:2*pi)*101,'k'); hold on;
scatter(tg2(:,1),tg2(:,2),cell2mat(tg2_c)',tg2_ts','Filled')
title('TgF344-AD Stops - Cue Absent')
axis image
axis off
colorbar
colormap(viridis(255))
subaxis(2,2,3)
plot(sin(0:2*pi/1000:2*pi)*101,cos(0:2*pi/1000:2*pi)*101,'k'); hold on;
scatter(wt1(:,1),wt1(:,2),cell2mat(wt1_c)',wt1_ts','Filled')
title('F344 Stops - Cue Present')
axis image
axis off
colorbar
colormap(viridis(255))
subaxis(2,2,4)
plot(sin(0:2*pi/1000:2*pi)*101,cos(0:2*pi/1000:2*pi)*101,'k'); hold on;
scatter(wt2(:,1),wt2(:,2),cell2mat(wt2_c)',wt2_ts','Filled')
title('F344 Stops - Cue Absent')
axis image
axis off
colorbar
colormap(viridis(255))
export_fig('D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Figures\PNGs\stops_duration_acrossTrial.png','-m4') 

%% Time spent moving for Movement segments 

tg1=horzcat(params.timeMoving{contains(param_idx,'Tg') & contains(param_idx,'D1')});
wt1=horzcat(params.timeMoving{contains(param_idx,'WT') & contains(param_idx,'D1')});
tg2=horzcat(params.timeMoving{contains(param_idx,'Tg') & contains(param_idx,'D2')});
wt2=horzcat(params.timeMoving{contains(param_idx,'WT') & contains(param_idx,'D2')});

fig=figure; 
fig.Color=[1 1 1];

AllStats  = stat_plot(cell2mat(wt1)',cell2mat(tg1)',{'WT','TG'},'Time Spent Moving (s)','plots',2)
title('Day 1 - Cue Present')
legend({'WT','Tg'})
AllStats  = stat_plot(cell2mat(wt2)',cell2mat(tg2)',{'WT','TG'},'Time Spent Moving (s)','plots',2)
title('Day 2 - Cue Absent')
legend({'WT','Tg'})

%% Path Lengths for Movement segments 

tg1=horzcat(params.segPL{contains(param_idx,'Tg') & contains(param_idx,'D1')});
wt1=horzcat(params.segPL{contains(param_idx,'WT') & contains(param_idx,'D1')});
tg2=horzcat(params.segPL{contains(param_idx,'Tg') & contains(param_idx,'D2')});
wt2=horzcat(params.segPL{contains(param_idx,'WT') & contains(param_idx,'D2')});

fig=figure; 
fig.Color=[1 1 1];

AllStats  = stat_plot(cell2mat(wt1)',cell2mat(tg1)',{'WT','TG'},'Length of Movment Segments (cm)','plots',2)
title('Day 1 - Cue Present')
legend({'WT','Tg'})
AllStats  = stat_plot(cell2mat(wt2)',cell2mat(tg2)',{'WT','TG'},'Length of Movment Segments (cm)','plots',2)
title('Day 2 - Cue Absent')
legend({'WT','Tg'})


%%
% Figure - Bar chart for quadrants - Path Length (cm)
tg1=vertcat(params.time_in_zone{contains(param_idx,'Tg') & contains(param_idx,'D1')});
wt1=vertcat(params.time_in_zone{contains(param_idx,'WT') & contains(param_idx,'D1')});
tg2=vertcat(params.time_in_zone{contains(param_idx,'Tg') & contains(param_idx,'D2')});
wt2=vertcat(params.time_in_zone{contains(param_idx,'WT') & contains(param_idx,'D2')});

stat_plot(wt2,tg2,{'WT','TG'},{'Time in Cue Zone'},'plots',2,'plottype','beeswarm')

subplot(1,2,1)
plotspread_wrapper(wt1,tg1,{'WT','TG'})
title('Time in Cue Zone - Cue Present')
ylabel('time (s)')
subplot(1,2,2)
plotspread_wrapper(wt2,tg2,{'WT','TG'})
title('Time in Cue Zone - Cue Absent')
ylabel('time (s)')

for i=1:size(params,1)
    fig=figure;
    fig.Color=[1 1 1];
    histogram(wt1,60)
    xlabel('Binned Degrees')
    title(param_idx{i})
end

%% Heading across entire trial 
tg1=vertcat(params.HD{contains(param_idx,'Tg') & contains(param_idx,'D1')});
wt1=vertcat(params.HD{contains(param_idx,'WT') & contains(param_idx,'D1')});
tg2=vertcat(params.HD{contains(param_idx,'Tg') & contains(param_idx,'D2')});
wt2=vertcat(params.HD{contains(param_idx,'WT') & contains(param_idx,'D2')});

fig=figure;
fig.Color=[1 1 1];
subplot(1,2,1)
h1=histogram(wt1,10);hold on
h2=histogram(tg1,10);
set(h1,'FaceColor',[.1 .1 .1],'EdgeColor','k')
set(h2,'FaceColor','r','EdgeColor','k')
set(gca,'box','off','FontWeight','bold','FontSize',18,'LineWidth',3)
xlabel('Binned Degrees')
subplot(1,2,2)
h1=histogram(wt2,30);hold on
h2=histogram(tg2,30);
set(h1,'FaceColor',[.1 .1 .1],'EdgeColor','k')
set(h2,'FaceColor','r','EdgeColor','k')
set(gca,'box','off','FontWeight','bold','FontSize',18,'LineWidth',3)
xlabel('Binned Degrees')


%% Heading across entire trial 
tg1=horzcat(params.timeMoving{contains(param_idx,'Tg') & contains(param_idx,'D1')});
wt1=horzcat(params.timeMoving{contains(param_idx,'WT') & contains(param_idx,'D1')});
tg2=horzcat(params.timeMoving{contains(param_idx,'Tg') & contains(param_idx,'D2')});
wt2=horzcat(params.timeMoving{contains(param_idx,'WT') & contains(param_idx,'D2')});

fig=figure;
fig.Color=[1 1 1];
subplot(1,2,1)
h1=histogram(cell2mat(wt1),30);hold on
h2=histogram(cell2mat(tg1),30);
set(h1,'FaceColor',[.1 .1 .1],'EdgeColor','k')
set(h2,'FaceColor','r','EdgeColor','k')
set(gca,'box','off','FontWeight','bold','FontSize',18,'LineWidth',3)
xlabel('Segments of Time Spent Moving')
subplot(1,2,2)
h1=histogram(cell2mat(wt2),30);hold on
h2=histogram(cell2mat(tg2),30);
set(h1,'FaceColor',[.1 .1 .1],'EdgeColor','k')
set(h2,'FaceColor','r','EdgeColor','k')
set(gca,'box','off','FontWeight','bold','FontSize',18,'LineWidth',3)
xlabel('Segments of Time Spent Moving')
    

AllStats  = stat_plot(cell2mat(wt1)',cell2mat(tg1)',{'WT','TG'},'Time of Movment Segments (s)','plots',2)
title('Day 1 - Cue Present')
legend({'WT','Tg'})
AllStats  = stat_plot(cell2mat(wt2)',cell2mat(tg2)',{'WT','TG'},'Time of Movment Segments (s)','plots',2)
title('Day 2 - Cue Absent')
legend({'WT','Tg'})


%% Heading across entire trial 
tg1=vertcat(params.overall_body_dir{contains(param_idx,'Tg') & contains(param_idx,'D1')});
wt1=vertcat(params.overall_body_dir{contains(param_idx,'WT') & contains(param_idx,'D1')});
tg2=vertcat(params.overall_body_dir{contains(param_idx,'Tg') & contains(param_idx,'D2')});
wt2=vertcat(params.overall_body_dir{contains(param_idx,'WT') & contains(param_idx,'D2')});

fig=figure;
fig.Color=[1 1 1];
subplot(1,2,1)
h1=histogram(wt1,30);hold on
h2=histogram(tg1,30);
set(h1,'FaceColor',[.1 .1 .1],'EdgeColor','k')
set(h2,'FaceColor','r','EdgeColor','k')
set(gca,'box','off','FontWeight','bold','FontSize',18,'LineWidth',3)
xlabel('Binned Degrees')
subplot(1,2,2)
h1=histogram(wt2,30);hold on
h2=histogram(tg2,30);
set(h1,'FaceColor',[.1 .1 .1],'EdgeColor','k')
set(h2,'FaceColor','r','EdgeColor','k')
set(gca,'box','off','FontWeight','bold','FontSize',18,'LineWidth',3)
xlabel('Binned Degrees')

