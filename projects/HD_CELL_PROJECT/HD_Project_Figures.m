% HD_Project_Figures
% this script works in combination with the HDdata structure output from
% HD_Project_TuningCurve_v2.m


%% HD CELL EXAMPLES

close all
fig=figure;fig.Color=[1 1 1];
scatter(HDdata.ATN.RLength,HDdata.ATN.within_Coeff,'filled');hold on
scatter(HDdata.PoSDeep.RLength,HDdata.PoSDeep.within_Coeff,'filled')
scatter(HDdata.PoSSup.RLength,HDdata.PoSSup.within_Coeff,'filled')
scatter(HDdata.MECDeepLayers.RLength,HDdata.MECDeepLayers.within_Coeff,'filled')
scatter(HDdata.MECSupLayers.RLength,HDdata.MECSupLayers.within_Coeff,'filled')
scatter(HDdata.PaSDeep.RLength,HDdata.PaSDeep.within_Coeff,'filled')
scatter(HDdata.PaSSup.RLength,HDdata.PaSSup.within_Coeff,'filled')
set(gca,'box','off','LineWidth',2,'FontWeight','bold','FontSize',20,'TickLength',[0;0])
xlabel('RLength')
ylabel('Stability (r)')
legend(areas,'Location','Best');

% print(fig,'-depsc','-tiff','-loose',['/Users/ryanharvey/Downloads',filesep,'RLengthbyStability_Fig.eps'])

close all
atnR=HDdata.ATN.RLength;
posR=[HDdata.PoSDeep.RLength;HDdata.PoSSup.RLength];
mecR=[HDdata.MECDeepLayers.RLength;HDdata.MECSupLayers.RLength];
pasR=[HDdata.PaSDeep.RLength;HDdata.PaSSup.RLength];

atnS=[HDdata.ATN.within_Coeff];
posS=[HDdata.PoSDeep.within_Coeff;HDdata.PoSSup.within_Coeff];
mecS=[HDdata.MECDeepLayers.within_Coeff;HDdata.MECSupLayers.within_Coeff];
pasS=[HDdata.PaSDeep.within_Coeff;HDdata.PaSSup.within_Coeff];


fig=figure; fig.Color=[1 1 1]; fig.OuterPosition=[1 6 960 1052];
findexamples(atnR,atnS,HDdata.ATN.rawTuningCurve);
print(fig,'-depsc','-tiff','-loose',['/Users/ryanharvey/Downloads',filesep,'atn_Fig.eps'])
close all

fig=figure; fig.Color=[1 1 1]; fig.OuterPosition=[1 6 960 1052];
findexamples(posR,posS,[HDdata.PoSDeep.rawTuningCurve;HDdata.PoSSup.rawTuningCurve]);
print(fig,'-depsc','-tiff','-loose',['/Users/ryanharvey/Downloads',filesep,'pos_Fig.eps'])
close all

fig=figure; fig.Color=[1 1 1]; fig.OuterPosition=[1 6 960 1052];
findexamples(mecR,mecS,[HDdata.MECDeepLayers.rawTuningCurve;HDdata.MECSupLayers.rawTuningCurve]);
print(fig,'-depsc','-tiff','-loose',['/Users/ryanharvey/Downloads',filesep,'mec_Fig.eps'])
close all

fig=figure; fig.Color=[1 1 1]; fig.OuterPosition=[1 6 960 1052];
findexamples(pasS,pasR,[HDdata.PaSDeep.rawTuningCurve;HDdata.PaSSup.rawTuningCurve]);
print(fig,'-depsc','-tiff','-loose',['/Users/ryanharvey/Downloads',filesep,'pas_Fig.eps'])
close all

%% BORDER SCORE
co=[ 0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840];
set(groot,'defaultAxesColorOrder',co)

close all
ATN=[HDdata.ATN.borderscore,HDdata.ATN.RLength];
PoS=[[HDdata.PoSDeep.borderscore;HDdata.PoSSup.borderscore],[HDdata.PoSDeep.RLength;HDdata.PoSSup.RLength]];
MEC=[[HDdata.MECDeepLayers.borderscore;HDdata.MECSupLayers.borderscore],[HDdata.MECDeepLayers.RLength;HDdata.MECSupLayers.RLength]];
PaS=[[HDdata.PaSDeep.borderscore;HDdata.PaSSup.borderscore],[HDdata.PaSDeep.RLength;HDdata.PaSSup.RLength]];

fig=figure;fig.Color=[1 1 1];
scatter(ATN(:,1),ATN(:,2),'filled');hold on
scatter(PoS(:,1),PoS(:,2),'filled')
scatter(MEC(:,1),MEC(:,2),'filled')
scatter(PaS(:,1),PaS(:,2),'filled')

set(gca,'box','off','LineWidth',2,'FontWeight','bold','FontSize',20,'TickLength',[0;0])
xlabel('Border Score')
ylabel('RLength')
legend({'ATN','PoS','MEC','PaS'},'Location','Best');
%% BORDER MODULATION
 co=[  0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840];
set(groot,'defaultAxesColorOrder',co)
close all;clear data
% data.ATN=max(HDdata.ATN.bordermodulation,[],2);
data.PoSDeep=max(HDdata.PoSDeep.bordermodulation,[],2);
data.PoSSup=max(HDdata.PoSSup.bordermodulation,[],2);
data.MECDeepLayers=max(HDdata.MECDeepLayers.bordermodulation,[],2);
data.MECSupLayers=max(HDdata.MECSupLayers.bordermodulation,[],2);
data.PaSDeep=max(HDdata.PaSDeep.bordermodulation,[],2);
data.PaSSup=max(HDdata.PaSSup.bordermodulation,[],2);

fig=figure;fig.Color=[1 1 1];
[freq]=ECDF_plot(data,'Border Modulation');
% print(fig,'-depsc','-tiff','-loose',['/Users/ryanharvey/Downloads',filesep,'bordermodulation.eps'])
print(fig,'-dmeta',['C:\Users\ryanh\Dropbox\school work\UNM\Lab\Projects\Multi_Region_HDcells',filesep,'bordermodulation.emf'])

%% BORDER SCORE TO BORDER MODULATION CORRELATION
PoSDeep=corr2(max(HDdata.PoSDeep.bordermodulation,[],2),HDdata.PoSDeep.borderscore)
figure;scatter(max(HDdata.PoSDeep.bordermodulation,[],2),HDdata.PoSDeep.borderscore)

PoSSup=corr2(max(HDdata.PoSSup.bordermodulation,[],2),HDdata.PoSSup.borderscore)
figure;scatter(max(HDdata.PoSSup.bordermodulation,[],2),HDdata.PoSSup.borderscore)

MECDeepLayers=corr2(max(HDdata.MECDeepLayers.bordermodulation,[],2),HDdata.MECDeepLayers.borderscore)
figure;scatter(max(HDdata.MECDeepLayers.bordermodulation,[],2),HDdata.MECDeepLayers.borderscore)

MECSupLayers=corr2(max(HDdata.MECSupLayers.bordermodulation,[],2),HDdata.MECSupLayers.borderscore)
figure;scatter(max(HDdata.MECSupLayers.bordermodulation,[],2),HDdata.MECSupLayers.borderscore)
 
PaSDeep=corr2(max(HDdata.PaSDeep.bordermodulation,[],2),HDdata.PaSDeep.borderscore)
figure;scatter(max(HDdata.PaSDeep.bordermodulation,[],2),HDdata.PaSDeep.borderscore)

PaSSup=corr2(max(HDdata.PaSSup.bordermodulation,[],2),HDdata.PaSSup.borderscore)
figure;scatter(max(HDdata.PaSSup.bordermodulation,[],2),HDdata.PaSSup.borderscore)

%% EGOCENTRIC MODULATION
% co=[ 0    0.4470    0.7410;
  co=[  0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840];
set(groot,'defaultAxesColorOrder',co)

close all;
clear data
% % data.ATN=max(HDdata.ATN.egocentricmodulation,[],2);
% data.ATN=nan(100);
data.PoSDeep=max(HDdata.PoSDeep.egocentricmodulation,[],2);
data.PoSSup=max(HDdata.PoSSup.egocentricmodulation,[],2);
data.MECDeepLayers=max(HDdata.MECDeepLayers.egocentricmodulation,[],2);
data.MECSupLayers=max(HDdata.MECSupLayers.egocentricmodulation,[],2);
data.PaSDeep=max(HDdata.PaSDeep.egocentricmodulation,[],2);
data.PaSSup=max(HDdata.PaSSup.egocentricmodulation,[],2);

fig=figure;fig.Color=[1 1 1];
[freq,p]=ECDF_plot(data,'Egocentric Modulation');

print(fig,'-dmeta',['C:\Users\ryanh\Dropbox\school work\UNM\Lab\Projects\Multi_Region_HDcells',filesep,'EgocentricModulation.emf'])
%% INFO CONTENT
close all;clear data
data.ATN=HDdata.ATN.informationContent;
data.PoSDeep=HDdata.PoSDeep.informationContent;
data.PoSSup=HDdata.PoSSup.informationContent;
data.MECDeepLayers=HDdata.MECDeepLayers.informationContent;
data.MECSupLayers=HDdata.MECSupLayers.informationContent;
data.PaSDeep=HDdata.PaSDeep.informationContent;
data.PaSSup=HDdata.PaSSup.informationContent;

fig=figure;fig.Color=[1 1 1];
[freq]=ECDF_plot(data,'Information Content (bits/spk)');

% print(fig,'-depsc','-tiff','-loose',['/Users/ryanharvey/Downloads',filesep,'bordermodulation.eps'])
print(fig,'-dmeta',['C:\Users\ryanh\Dropbox\school work\UNM\Lab\Projects\Multi_Region_HDcells',filesep,'InformationContent.emf'])

%% Grid Score
close all;clear data
data.ATN=HDdata.ATN.gridscore;
data.PoSDeep=HDdata.PoSDeep.gridscore;
data.PoSSup=HDdata.PoSSup.gridscore;
data.MECDeepLayers=HDdata.MECDeepLayers.gridscore;
data.MECSupLayers=HDdata.MECSupLayers.gridscore;
data.PaSDeep=HDdata.PaSDeep.gridscore;
data.PaSSup=HDdata.PaSSup.gridscore;

fig=figure;fig.Color=[1 1 1];
[freq]=ECDF_plot(data,'Grid Score');

% print(fig,'-depsc','-tiff','-loose',['/Users/ryanharvey/Downloads',filesep,'bordermodulation.eps'])
print(fig,'-dmeta',['C:\Users\ryanh\Dropbox\school work\UNM\Lab\Projects\Multi_Region_HDcells',filesep,'gridness.emf'])


%% THETA MODULATION FIGURES
clear data
data.ATN=HDdata.ATN.thetaindex;
data.PoSDeep=HDdata.PoSDeep.thetaindex;
data.PoSSup=HDdata.PoSSup.thetaindex;
data.MECDeepLayers=HDdata.MECDeepLayers.thetaindex;
data.MECSupLayers=HDdata.MECSupLayers.thetaindex;
data.PaSDeep=HDdata.PaSDeep.thetaindex;
data.PaSSup=HDdata.PaSSup.thetaindex;

fig=figure;fig.Color=[1 1 1];
[freq]=ECDF_plot(data,'Theta Index');
% print(fig,'-depsc', '-r600',['/Users/ryanharvey/Downloads',filesep,'ThetaModFiglayer.eps'])

%% THETA MODULATION FIGURES
clear data
data.ATN=HDdata.ATN.thetaindex;
data.PoS=[HDdata.PoSDeep.thetaindex;HDdata.PoSSup.thetaindex];
data.MEC=[HDdata.MECDeepLayers.thetaindex;HDdata.MECSupLayers.thetaindex];
data.PaS=[HDdata.PaSDeep.thetaindex;HDdata.PaSSup.thetaindex];

fig=figure;fig.Color=[1 1 1];
[freq]=ECDF_plot(data,'ThetaIndex');
print(fig,'-depsc', '-r600',['/Users/ryanharvey/Downloads',filesep,'ThetaModFigregion.eps'])

%% THETA MODULATED PROPORTIONS
% by layer
close all
y=[sum(HDdata.MECSupLayers.modulated==1)/length(HDdata.MECSupLayers.modulated),sum(HDdata.MECSupLayers.modulated==0)/length(HDdata.MECSupLayers.modulated);...
    sum(HDdata.PaSSup.modulated==1)/length(HDdata.PaSSup.modulated),sum(HDdata.PaSSup.modulated==0)/length(HDdata.PaSSup.modulated);
    sum(HDdata.MECDeepLayers.modulated==1)/length(HDdata.MECDeepLayers.modulated),sum(HDdata.MECDeepLayers.modulated==0)/length(HDdata.MECDeepLayers.modulated);...
    sum(HDdata.PaSDeep.modulated==1)/length(HDdata.PaSDeep.modulated),sum(HDdata.PaSDeep.modulated==0)/length(HDdata.PaSDeep.modulated);...
    sum(HDdata.PoSDeep.modulated==1)/length(HDdata.PoSDeep.modulated),sum(HDdata.PoSDeep.modulated==0)/length(HDdata.PoSDeep.modulated);...
    sum(HDdata.PoSSup.modulated==1)/length(HDdata.PoSSup.modulated),sum(HDdata.PoSSup.modulated==0)/length(HDdata.PoSSup.modulated);...
    sum(HDdata.ATN.modulated==1)/length(HDdata.ATN.modulated),sum(HDdata.ATN.modulated==0)/length(HDdata.ATN.modulated)];
categories={'MECSupLayers','PaSSup','MECDeepLayers','PaSDeep','PoSDeep','PoSSup','ATN'};

barfiglayer=figure;barfiglayer.Color=[1 1 1];
b=barh(y,'stacked')
xlabel('Proportion')
ylabel('Layer/Region')
set(b,'LineStyle','none');
set(gca,'yticklabel',categories,'box','off','LineWidth',2,'FontWeight','bold','FontSize',20,'TickLength',[0;0])
b(1, 1).FaceColor   = [0.929,  0.694,  0.125];
b(1, 2).FaceColor   = [0.3,  0.3,  0.3];

% by region
% close all
mec=[HDdata.MECSupLayers.modulated;HDdata.MECDeepLayers.modulated];
pas=[HDdata.PaSSup.modulated;HDdata.PaSDeep.modulated];
pos=[HDdata.PoSSup.modulated;HDdata.PoSDeep.modulated];
atn=[HDdata.ATN.modulated];

y=[sum(mec==1)/length(mec),sum(mec==0)/length(mec);...
    sum(pas==1)/length(pas),sum(pas==0)/length(pas);...
    sum(pos==1)/length(pos),sum(pos==0)/length(pos);...
    sum(atn==1)/length(atn),sum(atn==0)/length(atn)];
categories={'MEC','PaS','PoS','ATN'};

barfigregion=figure;barfigregion.Color=[1 1 1];
b=barh(y,'stacked')
xlabel('Proportion')
ylabel('Layer/Region')
set(b,'LineStyle','none');
set(gca,'yticklabel',categories,'box','off','LineWidth',2,'FontWeight','bold','FontSize',20,'TickLength',[0;0])
b(1, 1).FaceColor   = [0.929,  0.694,  0.125];
b(1, 2).FaceColor   = [0.3,  0.3,  0.3];

% print(barfiglayer,'-depsc','-tiff','-loose',['/Users/ryanharvey/Downloads',filesep,'ThetaModulationLayer_bar_Fig.eps'])

% print(barfigregion,'-depsc','-tiff','-loose',['/Users/ryanharvey/Downloads',filesep,'ThetaModulationRegion_bar_Fig.eps'])


%% POP VECTOR SORTED BY RLENGTH BY LAYERS
folders=fieldnames(HDdata);

for a=1:length(folders)
    tuning_Fig=figure; tuning_Fig.Color=[1 1 1];
    imagesc(HDdata.(folders{a}).TuningCurve);hold on
    %     plot(rescale(sum(HDdata.(folders{a}).TuningCurve),1,size(HDdata.(folders{a}).TuningCurve,1)),'w','LineWidth',3)
    axis xy;
    box off;
    ylabel('Cells')
    set(gca,'XTick',[1,30.7,60],'XTickLabel',[-180 0 180],'FontSize',20,'FontWeight','bold')
    title(folders{a})
    print(tuning_Fig,'-depsc', '-r600',['/Users/ryanharvey/Downloads',filesep,folders{a},'tuning_Fig.eps'])
    close all
end
%% DIAG POP VEC BY LAYERS

atn=HDdata.ATN.rawTuningCurve;
pos=[HDdata.PoSSup.rawTuningCurve;HDdata.PoSDeep.rawTuningCurve];
mec=[HDdata.MECSupLayers.rawTuningCurve;HDdata.MECDeepLayers.rawTuningCurve];
pas=[HDdata.PaSSup.rawTuningCurve;HDdata.PaSDeep.rawTuningCurve];

[atn]=arrangenorm(HDdata.ATN.rawTuningCurve);
[PoSDeep]=arrangenorm(HDdata.PoSDeep.rawTuningCurve);
[PoSSup]=arrangenorm(HDdata.PoSSup.rawTuningCurve);
[MECDeepLayers]=arrangenorm(HDdata.MECDeepLayers.rawTuningCurve);
[MECSupLayers]=arrangenorm(HDdata.MECSupLayers.rawTuningCurve);
[PaSDeep]=arrangenorm(HDdata.PaSDeep.rawTuningCurve);
[PaSSup]=arrangenorm(HDdata.PaSSup.rawTuningCurve);


data={atn,PoSDeep,PoSSup,MECDeepLayers,MECSupLayers,PaSDeep,PaSSup};
regions=fieldnames(HDdata);


for a=1:length(data)
    tuning_Fig=figure; tuning_Fig.Color=[1 1 1];
    imagesc(data{a});hold on
    axis xy;
    box off;
    ylabel('Cells')
    set(gca,'XTick',[1,30.7,60],'XTickLabel',[-180 0 180],'FontSize',20,'FontWeight','bold')
    title(regions{a})
    print(tuning_Fig,'-depsc', '-r600',['/Users/ryanharvey/Downloads',filesep,regions{a},'tuning_Fig.eps'])
    close all
end

%% POOL BRAIN REGION & CREATE POP VEC SORTED BY RLENGTH
atn=HDdata.ATN.rawTuningCurve;
pos=[HDdata.PoSSup.rawTuningCurve;HDdata.PoSDeep.rawTuningCurve];
mec=[HDdata.MECSupLayers.rawTuningCurve;HDdata.MECDeepLayers.rawTuningCurve];
pas=[HDdata.PaSSup.rawTuningCurve;HDdata.PaSDeep.rawTuningCurve];

[~,atnri]=sort(HDdata.ATN.RLength);
[~,posri]=sort([HDdata.PoSSup.RLength;HDdata.PoSDeep.RLength]);
[~,mecri]=sort([HDdata.MECSupLayers.RLength;HDdata.MECDeepLayers.RLength]);
[~,pasri]=sort([HDdata.PaSSup.RLength;HDdata.PaSDeep.RLength]);

atn=atn(atnri,:);
pos=pos(posri,:);
mec=mec(mecri,:);
pas=pas(pasri,:);

[atn]=shiftnorm(atn);
[pos]=shiftnorm(pos);
[mec]=shiftnorm(mec);
[pas]=shiftnorm(pas);

data={atn,pos,mec,pas};
regions={'atn','pos','mec','pas'};

for a=1:length(data)
    tuning_Fig=figure; tuning_Fig.Color=[1 1 1];
    imagesc(data{a});hold on
    axis xy;
    box off;
    ylabel('Cells')
    set(gca,'XTick',[1,30.7,60],'XTickLabel',[-180 0 180],'FontSize',20,'FontWeight','bold')
    title(regions{a})
    print(tuning_Fig,'-depsc', '-tiff','-loose',['/Users/ryanharvey/Downloads',filesep,regions{a},'tuning_Fig.eps'])
    close all
end
%% DIAG POP VECTORS BY REGION
atn=HDdata.ATN.rawTuningCurve;
pos=[HDdata.PoSSup.rawTuningCurve;HDdata.PoSDeep.rawTuningCurve];
mec=[HDdata.MECSupLayers.rawTuningCurve;HDdata.MECDeepLayers.rawTuningCurve];
pas=[HDdata.PaSSup.rawTuningCurve;HDdata.PaSDeep.rawTuningCurve];

[atn]=arrangenorm(atn);
[pos]=arrangenorm(pos);
[mec]=arrangenorm(mec);
[pas]=arrangenorm(pas);

data={atn,pos,mec,pas};
regions={'atn','pos','mec','pas'};

for a=1:length(data)
    tuning_Fig=figure; tuning_Fig.Color=[1 1 1];
    imagesc(data{a});hold on
    axis xy;
    box off;
    ylabel('Cells')
    set(gca,'XTick',[1,30.7,60],'XTickLabel',[-180 0 180],'FontSize',20,'FontWeight','bold')
    title(regions{a})
    print(tuning_Fig,'-depsc', '-tiff','-loose',['/Users/ryanharvey/Downloads',filesep,regions{a},'tuning_Fig.eps'])
    close all
end


%% LOCAL FUNCTIONS

% SHIFT PEAK TO CENTER
function [out]=shiftnorm(in)
[~,I]=max(in,[],2);
middlebin=round(median(1:size(in,2)));
for i=1:size(in,1)
    out(i,:)=rescale(circshift(in(i,:),(middlebin-I(i,:))-1,2),0,1);
end
end

% FOR DIAG POP VECTOR
function [matout]=arrangenorm(mat)
[~,I]=max(mat,[],2);
[~,I2]=sort(I);
matout=mat(I2,:);
for i=1:size(matout,1)
    matout(i,:)=rescale(matout(i,:),0,1);
end
end

function findexamples(atnR,atnS,tuning)
medR=median(atnR);
medS=median(atnS);
stdR=std(atnR);
stdS=std(atnS);

% TOP
topI=find(atnR>medR+stdR & atnS>medS+stdS);
topI=topI(randi([1,length(topI)],2,1));

% MIDDLE
middleI=find([atnR<medR+stdR & atnR>medR-stdR] & [atnS<medS+stdS & atnS>medS-stdS]);
middleI=middleI(randi([1,length(middleI)],2,1));

% BOTTOM
bottomI=find(atnR<medR-stdR & atnS<medS-stdS);
bottomI=bottomI(randi([1,length(bottomI)],2,1));

subplot(3,2,1)
plot(tuning(topI(1),:),'k','LineWidth',2)
ylabel('Firing Rate (hz)')
title(['RLength: ',num2str(atnR(topI(1))),' Stability: ',num2str(atnS(topI(1)))])
set(gca,'XTick',[1,30.7,60],'XTickLabel',[-180 0 180],'FontSize',20,'FontWeight','bold','TickLength',[0;0],'box','off','LineWidth',2)

subplot(3,2,2)
plot(tuning(topI(2),:),'k','LineWidth',2)
ylabel('Firing Rate (hz)')
title(['RLength: ',num2str(atnR(topI(2))),' Stability: ',num2str(atnS(topI(2)))])
set(gca,'XTick',[1,30.7,60],'XTickLabel',[-180 0 180],'FontSize',20,'FontWeight','bold','TickLength',[0;0],'box','off','LineWidth',2)

subplot(3,2,3)
plot(tuning(middleI(1),:),'k','LineWidth',2)
ylabel('Firing Rate (hz)')
title(['RLength: ',num2str(atnR(middleI(1))),' Stability: ',num2str(atnS(middleI(1)))])
set(gca,'XTick',[1,30.7,60],'XTickLabel',[-180 0 180],'FontSize',20,'FontWeight','bold','TickLength',[0;0],'box','off','LineWidth',2)

subplot(3,2,4)
plot(tuning(middleI(2),:),'k','LineWidth',2)
ylabel('Firing Rate (hz)')
title(['RLength: ',num2str(atnR(middleI(2))),' Stability: ',num2str(atnS(middleI(2)))])
set(gca,'XTick',[1,30.7,60],'XTickLabel',[-180 0 180],'FontSize',20,'FontWeight','bold','TickLength',[0;0],'box','off','LineWidth',2)

subplot(3,2,5)
plot(tuning(bottomI(1),:),'k','LineWidth',2)
ylabel('Firing Rate (hz)')
title(['RLength: ',num2str(atnR(bottomI(1))),' Stability: ',num2str(atnS(bottomI(1)))])
set(gca,'XTick',[1,30.7,60],'XTickLabel',[-180 0 180],'FontSize',20,'FontWeight','bold','TickLength',[0;0],'box','off','LineWidth',2)

subplot(3,2,6)
plot(tuning(bottomI(2),:),'k','LineWidth',2)
ylabel('Firing Rate (hz)')
title(['RLength: ',num2str(atnR(bottomI(2))),' Stability: ',num2str(atnS(bottomI(2)))])
set(gca,'XTick',[1,30.7,60],'XTickLabel',[-180 0 180],'FontSize',20,'FontWeight','bold','TickLength',[0;0],'box','off','LineWidth',2)

end