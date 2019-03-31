
%% Composite Ratemaps 

%Rotate heatmap relative to home base (primary) location 

% figure; subplot(1,2,1); plot(cos(0:2*pi/1000:2*pi)*101,sin(0:2*pi/1000:2*pi)*101,'-k');
% hold on; plot(params.transcoords{1}(:,1),params.transcoords{1}(:,2),'.k')
% hold on; plot(params.HBcenter{1}(:,1),params.HBcenter{1}(:,2),'*r')
pl=1;
pd=1;
sl=1;
sd=1;

for j=1:size(params.VideoName)
    
center=params.HBcenter{j};
XY=params.transcoords{j};

centerHB =[0; 0];
centerXY = repmat([0; 0], 1, length(XY));

%% compute degree of rotation (theta)
angle= atan2d(center(:,2)-0,center(:,1)-0) + 360*((center(:,2)-0)<0);
theta=deg2rad(360-angle);

%% define a counter-clockwise rotation matrix
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

% do the rotation...
rotHB = R*(center' - centerHB) + centerHB;% shift points in the plane so that the center of rotation is at the origin% apply the rotation about the origin
rotXY = R*([XY]' - centerXY) + centerXY;

rotXY=rotXY';
rotHB=rotHB';

% subplot(1,2,2); plot(cos(0:2*pi/1000:2*pi)*101,sin(0:2*pi/1000:2*pi)*101,'-k');
% hold on; plot(rotXY(:,1)',rotXY(:,2)','.k');
% hold on; plot(rotHB(:,1)',rotHB(:,2)','*r');

%recreate heatmap
xedge=linspace(-(params.dia(j)/2),(params.dia(j)/2),round(params.dia(j)/binsize));
yedge=linspace(-(params.dia(j)/2),(params.dia(j)/2),round(params.dia(j)/binsize));

[map] = histcounts2(rotHB(:,1)',rotHB(:,2)',xedge,yedge);
map=flipud(imrotate(map,90));
map=map/fr;
rotatedMap= imgaussfilt(map, 1.5);  %STORE THIS!

if params.group{j}==1 && contains(params.VideoName{j},'Light')
    PAEMap_L(:,:,pl)=rescale(rotatedMap);
    pl=pl+1;
elseif params.group{j}==1 && contains(params.VideoName{j},'dark')
    PAEMap_D(:,:,pd)=rescale(rotatedMap);
    pd=pd+1;
elseif params.group{j}==0 && contains(params.VideoName{j},'Light')
    SACMap_L(:,:,sl)=rescale(rotatedMap);
    sl=sl+1;
elseif params.group{j}==0 && contains(params.VideoName{j},'dark')
    SACMap_D(:,:,sd)=rescale(rotatedMap);
    sd=sd+1;
end


end
%% Mean occupancy maps 

fig=figure; fig.Color=[1 1 1];

subaxis(2,2,1, 'Spacing', 0.03, 'Padding', 0);
SAC1=PerfectCircRateMap(mean(SACMap_D,3),0);
imAlpha=ones(size(SAC1));
imAlpha(isnan(SAC1))=0;
imagesc(SAC1,'AlphaData',imAlpha);
title('Sacc Light')
axis xy; axis off; hold on; box off; axis image;
hold on

subaxis(2,2,2, 'Spacing', 0.03, 'Padding', 0);
SAC2=PerfectCircRateMap(mean(SACMap_L,3),0);
imAlpha=ones(size(SAC2));
imAlpha(isnan(SAC2))=0;
imagesc(SAC2,'AlphaData',imAlpha);
title('Sacc Dark')
axis xy; axis off; hold on; box off; axis image;

hold on
subaxis(2,2,3,'Spacing', 0.03, 'Padding', 0)
PAE1=PerfectCircRateMap(mean(PAEMap_D,3),0);
imAlpha=ones(size(PAE1));
imAlpha(isnan(PAE1))=0;
imagesc(PAE1,'AlphaData',imAlpha);
title('PAE Light')
axis xy; axis off; hold on; box off; axis image;
hold on

subaxis(2,2,4,'Spacing', 0.03, 'Padding', 0)
PAE2=PerfectCircRateMap(mean(PAEMap_L,3),0);
imAlpha=ones(size(PAE2));
imAlpha(isnan(PAE2))=0;
imagesc(PAE2,'AlphaData',imAlpha);
title('PAE Dark')
axis xy; axis off; hold on; box off; axis image;
hold on


%% Plot histograms normalized by peak for dwell time

dwellQuad_SAC1=vertcat(params.dwellQuad{cell2mat(params.group)==0 & contains(string(param_idx),'Light')});
dwellQuad_SAC2=vertcat(params.dwellQuad{cell2mat(params.group)==0 & contains(string(param_idx),'dark')});
dwellQuad_PAE1=vertcat(params.dwellQuad{cell2mat(params.group)==1 & contains(string(param_idx),'Light')});
dwellQuad_PAE2=vertcat(params.dwellQuad{cell2mat(params.group)==1 & contains(string(param_idx),'dark')});


%Find Peak bin and sort by peak 
[SAC1_peak,SAC1_idx]=max(dwellQuad_SAC1,[],2);
[SAC2_peak,SAC2_idx]=max(dwellQuad_SAC2,[],2);
[PAEl_peak,PAE1_idx]=max(dwellQuad_PAE1,[],2);
[PAE2_peak,PAE2_idx]=max(dwellQuad_PAE2,[],2);

 middlebin=round(median(1:length(dwellQuad_SAC1)));
 for a=1:size(dwellQuad_SAC1,1)
 SAC1_shift(a,:)=circshift(dwellQuad_SAC1(a,:),(middlebin-SAC1_idx(a,1))-1);
 end

 for a=1:size(dwellQuad_SAC1,1)
 SAC2_shift(a,:)=circshift(dwellQuad_SAC2(a,:),(middlebin-SAC2_idx(a,1))-1);
 end

 for a=1:size(dwellQuad_SAC1,1)
 PAE1_shift(a,:)=circshift(dwellQuad_PAE1(a,:),(middlebin-PAE1_idx(a,1))-1);
 end

 for a=1:size(dwellQuad_SAC1,1)
 PAE2_shift(a,:)=circshift(dwellQuad_PAE2(a,:),(middlebin-PAE2_idx(a,1))-1);
 end

%Histograms for quadrant dwell time
f=figure; f.Color=[1 1 1];
subplot(1,2,1)
bar(sum(SAC1_shift,1)/sum(SAC1_shift(:)),'FaceColor','k','FaceAlpha',.5);
title('Quadrant Dwell Time for PAE (red) and SACC (black) in Dark')
ylim([0 .4])
xlabel('Quadrant Angle Centered at peak Bin')
ylabel('Normalized Dwell Time')
 hold on 
 bar(sum(PAE1_shift,1)/sum(PAE1_shift(:)),'FaceColor','r','FaceAlpha',.5);
 set(gca,'FontSize',14,'FontWeight','bold','FontName','Calibri')
box off 

subplot(1,2,2)
bar(sum(SAC2_shift,1)/sum(SAC2_shift(:)),'FaceColor','k','FaceAlpha',.5);
title('Quadrant Dwell Time for PAE (red) and SACC (black) in Light')
ylim([0 .4])
xlabel('Quadrant Angle Centered at peak Bin')
ylabel('Normalized Dwell Time')
 hold on 
 bar(sum(PAE2_shift,1)/sum(PAE2_shift(:)),'FaceColor','r','FaceAlpha',.5);
  set(gca,'FontSize',14,'FontWeight','bold','FontName','Calibri')
box off 

