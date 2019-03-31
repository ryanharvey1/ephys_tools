% LandMarkControl
clear;close all;clc
addpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis')

load('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/Ryan_et_al_2018_MasterData.mat')

controldisp=controlplacecelldata(:,14:15,2);
tiltdisp=tiltedplacecelldata(:,14:15,2);

load('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/Figures_For_Paper/Displacement_Shuffled_correlations.mat')

controlplacecelldata(controldisp(:,2)<prctile(rSHUFF,95),:,:)=[];
tiltedplacecelldata(tiltdisp(:,2)<prctile(rSHUFF,95),:,:)=[];


controldisp(controldisp(:,2)<prctile(rSHUFF,95),:)=[];
tiltdisp(tiltdisp(:,2)<prctile(rSHUFF,95),:)=[];


%% PIE CHART FOR AMOUNT OF CELLS THAT ROTATE
% don't rotate
norotcon=(sum(controldisp(:,1)>360-45 | controldisp(:,1)<45)/length(controldisp))*100
norottilt=(sum(tiltdisp(:,1)>360-45 | tiltdisp(:,1)<45)/length(tiltdisp))*100
% rotate with cue
rotcon=(sum(controldisp(:,1)>90-45 & controldisp(:,1)<90+45)/length(controldisp))*100
rottilt=(sum(tiltdisp(:,1)>90-45 & tiltdisp(:,1)<90+45)/length(tiltdisp))*100
% rotates not with cue 
otherrotcon=(sum(controldisp(:,1)>135 & controldisp(:,1)<315)/length(controldisp))*100
otherrottilt=(sum(tiltdisp(:,1)>135 & tiltdisp(:,1)<315)/length(tiltdisp))*100


pie_displacement_fig=figure; pie_displacement_fig.Color=[1 1 1];
ax=subplot(1,2,1);p1=pie([norotcon,rotcon,otherrotcon]);
title(ax,'Control');ax.FontSize=20;
hold on
ax=subplot(1,2,2);p2=pie([norottilt,rottilt,otherrottilt]);
title(ax,'Tilted');ax.FontSize=20;
%% CALCULATE SPATIAL MEASURES ON CELLS THAT DID AND DID NOT ROTATE WITH CUE
controlplacecellpaths(controldisp(:,1)>360-45 | controldisp(:,1)<45 | controldisp(:,1)>135 & controldisp(:,1)<315,:,:)
tiltedplacecellpaths(tiltdisp(:,1)>360-45 | tiltdisp(:,1)<45 | tiltdisp(:,1)>135 & tiltdisp(:,1)<315,:,:)

% don't rotate
norotcon=controlplacecelldata(controldisp(:,1)>360-45 | controldisp(:,1)<45 | controldisp(:,1)>135 & controldisp(:,1)<315,:,:);
norottilt=tiltedplacecelldata(tiltdisp(:,1)>360-45 | tiltdisp(:,1)<45 | tiltdisp(:,1)>135 & tiltdisp(:,1)<315,:,:);

% rotate with cue
rotcon=controlplacecelldata(controldisp(:,1)>90-45 & controldisp(:,1)<90+45,:,:);
rottilt=tiltedplacecelldata(tiltdisp(:,1)>90-45 & tiltdisp(:,1)<90+45,:,:);

[ controldifferences ] = ScatterBox(rotcon(:,:,2),norotcon(:,:,2),{'Control Rotate', 'Control No Rotate'},Varnames,1)

[ tilteddifferences ] = ScatterBox(rottilt(:,:,2),norottilt(:,:,2),{'Tilted Rotate', 'Tilted No Rotate'},Varnames,1)


