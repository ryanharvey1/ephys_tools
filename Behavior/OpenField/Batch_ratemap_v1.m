cd ('d:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\PAE_OF\TrackData\')
file=struct2table(dir( '**/*.txt'));
t=table2cell(file(:,1:2));
fps=30;
params.dia=202;
binsize=6;

for i=1:length(t)
    
    delimiterIn = '\t';
    headerlinesIn = 1;
    
    coords1 = importdata(t{i,1}, delimiterIn, headerlinesIn);
    
    tempX=coords1.data(:,1);
    tempY=coords1.data(:,2);
    data.ts{i}=linspace(0,size(tempX,1)/fps, size(tempX,1))';
    
    [procX,procY]=FixPos(tempX,tempY,data.ts{i});
    
    
    
    data.id{i}=t{i,1};
    data.x{i}=procX;
    data.y{i}=procY;
    
    clear coords1
    
end

%GET MAX/MIN
tempX=[];
tempY=[];
for ii=1:length(data.x)
    tempX=[tempX;data.x{ii}];
    tempY=[tempY;data.y{ii}];
end

%Get Params
xMax=max(tempX); xMin=min(tempX); yMax=max(tempY); yMin=min(tempY);

xedge=linspace(xMin,xMax,round(mazesize/binsize));
yedge=linspace(yMin,yMax,round(mazesize/binsize));

%Get BINS
for ii=1:length(data.x)
    % set up figure
    fig=figure;fig.Color=[1 1 1];
    
    subplot(8,3,1)
    plot(data.x{ii},data.y{ii},'.k')
    axis image
    axis off
    title(sprintf('session dur %4.2f (min)',data.ts{ii}(end)/60))
    
    subplot(8,3,2)
    [map] = histcounts2(data.x{ii},data.y{ii},xedge,yedge);
    map=flipud(imrotate(map,90));
    map=map/fps;
    
    data.map{ii}=map;
    
    map2=map;
    map2((map2)==0)=NaN;
    ax = gca;
    imAlpha=ones(size(map2));
    imAlpha(isnan(map2))=0;
    imagesc(map2,'AlphaData',imAlpha);
    axis xy; axis off; hold on; box off; axis image;
    colormap(ax,jet(255))
    title(sprintf('max occ %4.2f (min)',max(map(:))/60))

    %Plot Ratemap
    subplot(8,3,3)
    map3=PerfectCircRateMap(map,0);
    ax = gca;
    imAlpha=ones(size(map3));
    imAlpha(isnan(map3))=0;
    imagesc(map3,'AlphaData',imAlpha);
    axis xy; axis off; hold on; box off; axis image;
    colormap(ax,jet(255))
    title(sprintf('max occ %4.2f (min)',max(map(:))/60))
    
    jumps=0:300:1800;
    
    figIt=4;
    for iii=1:length(jumps)
        
        if iii==7
            continue
        end
        
        idx=data.ts{ii}>jumps(iii) & data.ts{ii}<jumps(iii+1);
        tempx=data.x{ii}(idx);
        tempy=data.y{ii}(idx);
        
        subplot(8,3,figIt)
        plot(tempx,tempy,'.k')
        axis image
        axis off
        title(sprintf('session dur %4.2f (min)',(length(tempx)/fps)/60))
        
        subplot(8,3,figIt+1)
        [map] = histcounts2(tempx,tempy,xedge,yedge);
        map=flipud(imrotate(map,90));
        map=map/fps;
        
        map2=map;
        map2((map2)==0)=NaN;
        ax = gca;
        imAlpha=ones(size(map2));
        imAlpha(isnan(map2))=0;
        imagesc(map2,'AlphaData',imAlpha);
        axis xy; axis off; hold on; box off; axis image;
        colormap(ax,jet(255))
        title(sprintf('max occ %4.2f (min)',max(map(:))/60))
        
        subplot(8,3,figIt+2)
        map3=PerfectCircRateMap(map,0);
        ax = gca;
        imAlpha=ones(size(map3));
        imAlpha(isnan(map3))=0;
        imagesc(map3,'AlphaData',imAlpha);
        axis xy; axis off; hold on; box off; axis image;
        colormap(ax,jet(255))
        title(sprintf('max occ %4.2f (min)',max(map(:))/60))
        
        figIt=figIt+3;
    end
    
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    print(gcf,'-dpng', '-r300',...
        ['C:\Users\Ben Clark''s Lab\Desktop\FASD\track_data_figures\',erase(data.id{ii},'.txt'),'.png'])
    close all
end

