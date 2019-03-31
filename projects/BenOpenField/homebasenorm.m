% homebasenorm
clear;close all;clc
S=importdata('/Users/RyanHarvey/Downloads/HomeBaseStudy_SmoothedScaled.xlsx');
getParam=readtable('/Users/RyanHarvey/Downloads/OpenFieldParamsSmoothedScaled.xlsx');
rats=fieldnames(S);
getParam(:,1)=[];
getParam=table2array(getParam);

T1=zeros(30,30);
T2=zeros(30,30);
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
    
    angle=atan2d(homebase(2)-8.5,homebase(1)-8.5) + 360*((homebase(2)-8.5)<0);
    
    rotateangle=(360-angle)-180;
    
    Omatrix=imrotate(Omatrix,-rotateangle);
    
    normA=Omatrix-min(Omatrix(:));
    Omatrix=normA./max(normA(:));
    
    Omatrix=imresize(Omatrix,(30/length(Omatrix)));
    
    if mod(j,2)==1
        T1=sum(cat(3,T1,Omatrix),3);
    elseif mod(j,2)==0
        T2=sum(cat(3,T2,Omatrix),3);
    end
end
subplot(1,2,1)
upsamRateMap=PerfectCircRateMap(T1,0);
imAlpha=ones(size(upsamRateMap));
imAlpha(isnan(upsamRateMap))=0;
figure(1)
imagesc(upsamRateMap,'AlphaData',imAlpha);hold on
set(gca,'color',[1 1 1]);
box off
axis image
axis off
colormap jet
shading interp
set(figure(1),'Position',[1 5 720 800]);

upsamRateMap=PerfectCircRateMap(T2,0);
imAlpha=ones(size(upsamRateMap));
imAlpha(isnan(upsamRateMap))=0;
subplot(1,2,2)
imagesc(upsamRateMap,'AlphaData',imAlpha);
box off
axis image
axis off
colormap jet
shading interp
print(figure(1),'-bestfit', '-dpdf', '-r600',char(strcat('/Users/RyanHarvey/Downloads',filesep,'BEN_HOMEBASE.pdf')))

