% Tilted_waveformproperties

% filter down to just principle cells
find(~contains(allcontrolpaths(:,:,1),controlplacecellpaths))
% 
% 
% * NEVER FINISHED *
% 
% clear;close all;clc
controldata(controldata(:,4,1)>10,:);
allcontrolpaths(controldata(:,4,1)>10 & controldata(:,5,1)<200,:);
%%
close all



InterneuronFig=figure;InterneuronFig.Color=[1 1 1];
p1=plot(wave(1,:),wave(6,:)*1000,'k');hold on

[v,vI]=min(wave(6,:));
[p,pI]=max(wave(6,:));
dur=plot([wave(1,pI),wave(1,vI)],[-90,-90],'k')
set(dur,'LineWidth',3)
axis off
% ylabel('?V')
% xlabel('ms')
title('Interneuron')
set(p1,'LineWidth',3)
set(gca,'box','off','LineWidth',2,'FontSize',20,'FontWeight','bold')
print(InterneuronFig,'-dpng', '-r600',['/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/Figures_For_Paper/InterneuronFig.png'])



PrincipleFig=figure;PrincipleFig.Color=[1 1 1];
p2=plot(principlewave(1,:),principlewave(3,:)*1000,'k');hold on

[v,vI]=min(principlewave(3,:));
[p,pI]=max(principlewave(3,:));
dur=plot([principlewave(1,pI),principlewave(1,vI)],[-90,-90],'k')
set(dur,'LineWidth',3)
axis off

% ylabel('?V')
% xlabel('ms')
title('Pyramidal')
set(p2,'LineWidth',3)
set(gca,'box','off','LineWidth',2,'FontSize',20,'FontWeight','bold')
    
print(PrincipleFig,'-dpng', '-r600',['/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/Figures_For_Paper/PrincipleFig.png'])


% propfig=figure;propfig.Color=[1 1 1]
% p3=scatter([controldata(:,4,1);tilteddata(:,4,1)],[controldata(:,5,1);tilteddata(:,5,1)],'filled')
% xlabel('Average Firing Rate (hz)')
% ylabel('Spike Width')
% set(p3,'MarkerFaceColor','k')


%%
controldata(controldata(:,4,1)<10,:);
allcontrolpaths(controldata(:,4,1)<10 & controldata(:,5,1)>250,:);



% LOAD SS3D TT FILE
% fn='/Users/RyanHarvey/Desktop/TT4 - MARGINAL2_2.ntt';
addpath(genpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/BClarkToolbox/Analyses/spikeCode/MEX'));
fn= '/Volumes/Ryan_4TB/Place_Cell_Data/Place Cells - Tilted Mice/2013-08-01_10-31-31 - SR011 - 1 Place cell 1Placey cell RY analyzed/TT1 - RY.ntt';
[Timestamps, ScNumbers, CellNumbers, Features, Samples, Header] =Nlx2MatSpike_v3(fn,[1 1 1 1 1],1,1,[]);
cell1=Samples(:,:,CellNumbers==1);
cell2=Samples(:,:,CellNumbers==2);
cell3=Samples(:,:,CellNumbers==3);
cell4=Samples(:,:,CellNumbers==4);
cell5=Samples(:,:,CellNumbers==5);
cell6=Samples(:,:,CellNumbers==6);
cell7=Samples(:,:,CellNumbers==7);
cell8=Samples(:,:,CellNumbers==8);

[r,c,d]=size(cell1);
waveforms=figure; waveforms.Color=[1 1 1];
for ii=1:4
    subplot(2,2,ii)
    for i=1:round(d/4)
        plot(cell1(:,ii,i))
        hold on
    end
end

% average wave forms
avgwave=[mean(cell1(:,1,:),3)';mean(cell1(:,2,:),3)';mean(cell1(:,3,:),3)';mean(cell1(:,4,:),3)'];
