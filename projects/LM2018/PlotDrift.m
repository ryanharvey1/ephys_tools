
% %CELL 31
% N=data.LB04.S20170421141652.drift.binnedFR{15, 1};  
% meanHD=data.LB04.S20170421141652.drift.meanHD{15, 1} ;
% plotIdx=data.LB04.S20170421141652.drift.plotIdx{15, 1};  
% 
% %CELL 2
% N=data.LB02.S2016082519444021.drift.binnedFR{2, 1};   
% meanHD=data.LB02.S2016082519444021.drift.meanHD{2, 1};  
% plotIdx=data.LB02.S2016082519444021.drift.plotIdx{2, 1}; 

for hdcell=1:length(TGids)
    N=data.(TGids{hdcell,1}).(TGids{hdcell,2}).drift.binnedFR{TGids{hdcell,3},1};
    meanHD=data.(TGids{hdcell,1}).(TGids{hdcell,2}).drift.meanHD{TGids{hdcell,3}, 1};
    plotIdx=data.(TGids{hdcell,1}).(TGids{hdcell,2}).drift.plotIdx{TGids{hdcell,3}, 1};
    fig1=figure; fig1.Color=[1 1 1]; plot(meanHD,'-','Color',[.75 .75 .75],'LineWidth',3); hold on;
    y=1:length(meanHD);
    scatter(y(plotIdx==1),meanHD(plotIdx==1,:),50,'r','o','filled');
    box off
    set(gca,'FontWeight','Bold','FontSize',20,'LineWidth',3)
    print(figure(fig1), '-dpng', '-r600',['d:\Users\BClarkLab\Desktop\Laura Temp\SpikebyAngle',filesep,['TG_spikeByAngle_', num2str(hdcell),'.png']])
close all
end

%% Tuning Curves

binSize=6;
Npoints     = round(10*6/binSize); %for 6degree bins
gw          = gausswin(Npoints,5); %alpha = 0.05
gw          = gw/sum(gw);

%Plot Tuning Curves over quarters 
for hdcell=1:length(TGids)
    fig=figure; fig.Color=[1 1 1];
    for quarter=1:4
        tempTun=data.(TGids{hdcell,1}).(TGids{hdcell,2}).InterTrialStability{TGids{hdcell,3}, 1}.hdTuning{quarter, 1};
        %smooth curve
        l           = length(tempTun);
        hdTmp       = [tempTun tempTun tempTun];
        hdTmp       = convn(hdTmp(:),gw,'same');
        hdSmoothed  = hdTmp(l+1:2*l);
        h=plot(hdSmoothed); hold on;
        if quarter==1
            h.LineWidth=3;
            h.Color='k';
        elseif quarter ==2
            h.LineWidth=3;
            h.Color='b';
        elseif quarter ==3
            h.LineWidth=3;
            h.Color='r';
        elseif quarter==4
            h.LineWidth=3;
            h.Color=[.75 .75 .75];
        end
        box off
        set(gca,'FontWeight','Bold','FontSize',20,'LineWidth',3)
    end
            print(figure(fig), '-dpng', '-r600',['d:\Users\BClarkLab\Desktop\Laura Temp\WithinStability',filesep,['TG_within_', num2str(hdcell),'.png']])
            close all
end