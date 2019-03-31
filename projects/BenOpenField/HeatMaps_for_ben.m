% HeatMaps_for_ben
clear;close all;clc
S=importdata('/Users/RyanHarvey/Downloads/HomeBaseStudy_SmoothedScaled.xlsx');
getParam=readtable('/Users/RyanHarvey/Downloads/OpenFieldParamsSmoothedScaled.xlsx'); 
rats=fieldnames(S);
getParam(:,1)=[];
getParam=table2array(getParam);

for j=1:length(getParam)
    XY=S.(rats{j});
    nBinsx = round(122/7); nBinsy = round(122/7);
    x_max=getParam(j,2);x_min=getParam(j,1);y_max=getParam(j,4);y_min=getParam(j,3);
    edges{1} = linspace(y_min, y_max, nBinsy+1);
    edges{2} = linspace(x_min, x_max, nBinsx+1);
    
    Omatrix = hist3([XY(:,2) XY(:,1)],'Edges',edges);
    Omatrix(end,:) = [];
    Omatrix(:,end) = [];
    
    % scale to size of binned data
    XYtemp=[XY;[x_max,y_max;x_min,y_min];[getParam(j,5),getParam(j,6)]];
    XYtemp(:,1)=rescale(XYtemp(:,1),1,length(Omatrix));
    XYtemp(:,2)=rescale(XYtemp(:,2),1,length(Omatrix));
    homebase=XYtemp(end,:);
    XYtemp([end-2:end],:)=[];
    XY=XYtemp;
    
    [ymaxrate,xmaxrate]=find(Omatrix==max(max(Omatrix)));
    d(j)=sqrt((xmaxrate-homebase(1))^2+(ymaxrate-homebase(2))^2);
    
    figure(j)
    Omatrix(Omatrix==0)=NaN;
    imAlpha=ones(size(Omatrix));
    imAlpha(isnan(Omatrix))=0;
    imagesc(Omatrix,'AlphaData',imAlpha); hold on;
    set(gca,'color',[1 1 1]);
    box off
    axis image
    axis on
    colormap jet
    shading interp
    set(figure(j),'Position',[1 5 720 800]);
    plot([xmaxrate homebase(1)], [ymaxrate homebase(2)],'Color','k','LineWidth',6); % to plot min line using indexed min
    print(figure(j),'-dpng','-r600',char(strcat('/Users/RyanHarvey/Downloads',filesep,rats{j},'.png')))
    close all
    
    
%     upsamRateMap=PerfectCircRateMap(Omatrix,0);
% %     [ymaxrateup,xmaxrateup]=find(upsamRateMap==max(max(upsamRateMap)));
% %     homebaseup=[homebase(1)*(xmaxrateup/xmaxrate),homebase(2)*(ymaxrateup/ymaxrate)];
%     imAlpha=ones(size(upsamRateMap));
%     imAlpha(isnan(upsamRateMap))=0;
%     figure(j)
%     imagesc(upsamRateMap,'AlphaData',imAlpha);
%     set(gca,'color',[1 1 1]);
%     box off
%     axis image
%     axis off
%     colormap jet
%     shading interp
%     set(figure(j),'Position',[1 5 720 800]);
%     print(figure(j),'-bestfit', '-dpdf', '-r600',char(strcat('/Users/RyanHarvey/Downloads',filesep,rats{j},'.pdf')))
%     close all
end
    d=d'*7;
