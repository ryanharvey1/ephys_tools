% TiltedOccupancy
% to show that tilted and control mice had similar occupancy
clear;clc;close all
load('/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/NewOCC_Map_workspace2.mat')

FigureLocation='/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/Figures_For_Paper';

occC=[];
sessions=fieldnames(OCCC);
for i=1:length(sessions)
    occC=[occC;OCCC.(sessions{i})];
end

occT=[];
sessions=fieldnames(OCCT);
for i=1:length(sessions)
    occT=[occT;OCCT.(sessions{i})];
end


%% PROGRESSIVELY MOVE IN AND CALCULATE % OF POINTS OUTSIDE POLYGON 
occMatrixC=occC;
occMatrixT=occT;
xmax=max(occMatrixC(:,1));
xmin=min(occMatrixC(:,1));
ymax=max(occMatrixC(:,2));
ymin=min(occMatrixC(:,2));
xC=median([xmax,xmin]); yC=median([ymax,ymin]); radC=(max(xmax-xmin,ymax-ymin))/2;

xmax=max(occMatrixT(:,1));
xmin=min(occMatrixT(:,1));
ymax=max(occMatrixT(:,2));
ymin=min(occMatrixT(:,2));

xT=median([xmax,xmin]); yT=median([ymax,ymin]); radT=(max(xmax-xmin,ymax-ymin))/2;
th = 0:pi/179.5:2*pi; % 0 to 2*pi(6.28318530717959) at 0.0175 increments to equal 360 points
% xunit = rad * cos(th) + x; yunit = rad * sin(th) + y;
ii=1;
clear controlper tiltper
close all
clc
% fig1=figure(1);plot(occC(:,1),occC(:,2),'LineWidth',1,'color','k');fig1.Color=[1 1 1];axis off;axis image;box off;hold on
% fig2=figure(2);plot(occT(:,1),occT(:,2),'LineWidth',1,'color','k');fig2.Color=[1 1 1];axis off;axis image;box off;hold on
% pause(.00001)
cirC=0:8:round(radC);
cirT=0:8:round(radT);
for i=1:18
    disp(['progress: ',num2str(i)])
    temprad=round(radC)-cirC(i);
    xunit=temprad*cos(th)+xC; 
    yunit=temprad*sin(th)+yC;
%     figure(1); plot(xunit,yunit,'LineWidth',1,'color','r'); hold on;   pause(.00001)
    [in,~]=inpolygon(occMatrixC(:,1),occMatrixC(:,2),xunit,yunit);
    controlper(ii)=(sum(~in)/length(occMatrixC))*100;
    
    temprad=round(radT)-cirT(i);
    xunit=temprad*cos(th)+xT; 
    yunit=temprad*sin(th)+yT;
%     figure(2); plot(xunit,yunit,'LineWidth',1,'color','r'); hold on;   pause(.00001)
    [in,~]=inpolygon(occMatrixT(:,1),occMatrixT(:,2),xunit,yunit);
    tiltper(ii)=(sum(~in)/length(occMatrixT))*100;
    ii=ii+1;
end
%% plot results from above
fig=figure;fig.Color=[1 1 1];

p1=plot(controlper);hold on
p2=plot(tiltper);

set(p1,'LineWidth',4)
set(p2,'LineWidth',4)

box off
legend({'Control','Tilted'},'FontSize',12,'Location','best','box','off','FontSize',15,'FontWeight','Bold')
xlabel('Distance From Wall');ylabel('Percent Occupancy')
ax=gca;
set(ax,'FontSize',20,'FontWeight','bold','LineWidth',2)

print(figure(5),'-bestfit','-dpdf', '-r600',[FigureLocation,filesep,'cummOccupancy.pdf'])

[ AllStats ] = ScatterBox(controlper',tiltper',{'Control','Tilted'},{'Percent_Occupancy'},2)
[p,h,stats] = ranksum(controlper',tiltper')
[ AllStats ] = CDFplots(controlper',tiltper',{'Control','Tilted'},{'Percent_Occupancy'},2)


%%
occMatrixC=occC;
% 
% occMatrixC(occMatrixC(:,1)>median(occC(:,1))+std(occC(:,1))*2,:)=[];
% occMatrixC(occMatrixC(:,1)<median(occC(:,1))-std(occC(:,1))*2,:)=[];
% 
% occMatrixC(occMatrixC(:,2)>median(occC(:,2))+std(occC(:,2))*2,:)=[];
% occMatrixC(occMatrixC(:,2)<median(occC(:,2))-std(occC(:,2))*2,:)=[];
% 
% 
    nBinsx = round(61/2.44); nBinsy = round(61/2.44);
% %     mat4occ=unique(occMatrix,'rows');
    MinY = min(occMatrixC(:,2));
    MaxY = max(occMatrixC(:,2));
    MinX = min(occMatrixC(:,1));
    MaxX = max(occMatrixC(:,1));
    edges{1} = linspace(MinY, MaxY, nBinsy+1);
    edges{2} = linspace(MinX, MaxX, nBinsx+1);
    
    Omatrix = hist3([occMatrixC(:,2) occMatrixC(:,1)],'Edges',edges);
    Omatrix(end,:)=[];
    Omatrix(:,end)=[];
    
    x=histcounts2(occMatrixC(:,2), occMatrixC(:,1),[nBinsx nBinsy]);
    occ=x./60;
    
    max(max(occ))
    min(min(occ(occ>0)))

    occ = Omatrix./60;
%     
    PerfectCircRateMap(Omatrix,1);
%      print(figure(1),'-dpng', '-r600',[FigureLocation,filesep,'occmapControl.png'])
% 
%              autocorr=SpatialAutoCorr(Omatrix,length(Omatrix));

%%
occMatrixC=occT;    
% 
% occMatrixC(occMatrixC(:,1)>median(occC(:,1))+std(occC(:,1))*2,:)=[];
% occMatrixC(occMatrixC(:,1)<median(occC(:,1))-std(occC(:,1))*2,:)=[];
% 
% occMatrixC(occMatrixC(:,2)>median(occC(:,2))+std(occC(:,2))*2,:)=[];
% occMatrixC(occMatrixC(:,2)<median(occC(:,2))-std(occC(:,2))*2,:)=[];
% % 
    nBinsx = round(61/2.44); nBinsy = round(61/2.44);
%     mat4occ=unique(occMatrix,'rows');
    MinY = min(occMatrixC(:,2));
    MaxY = max(occMatrixC(:,2));
    MinX = min(occMatrixC(:,1));
    MaxX = max(occMatrixC(:,1));
    edges{1} = linspace(MinY, MaxY, nBinsy+1);
    edges{2} = linspace(MinX, MaxX, nBinsx+1);
%     
    OmatrixT = hist3([occMatrixC(:,2) occMatrixC(:,1)],'Edges',edges);
     OmatrixT(end,:)=[];
    OmatrixT(:,end)=[];
% 
% 
%     occt = OmatrixT/60;
%     
    PerfectCircRateMap(OmatrixT,1);
%          print(figure(2),'-dpng', '-r600',[FigureLocation,filesep,'occmapTilted.png'])
% 
%     
