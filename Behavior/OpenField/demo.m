
%     subtraction = (params.dia(i)/2)/3; %%calculates a third of the radius to use to plot the snaller circle zones in line__
    x0=0;%%origin
    y0=0;%%origin
    r=params.dia(34)/2;%%radius
    teta=-pi:0.01:pi;
    x=r*cos(teta)+x0;
    y=r*sin(teta)+y0;
    figure
    plot(x,y)
    hold on
    scatter(x0,y0,'or')
    axis square 
%     ----------------------------------------
%     divide your circle to n sectors
    n=16;
    tet=linspace(-pi,pi,n+1); %%allows for even space between each line??
    xi=(r+r*.25)*cos(tet)+x0;
    yi=(r+r*.25)*sin(tet)+y0;
    
    hold on;
    
%         for k=1:numel(xi)
        plot([x0 xi(16)],[y0 yi(16)])
        hold on
        plot([x0 xi(17)],[y0 yi(17)])
        hold on
        plot(params.transcoords{34,1}(:,1),params.transcoords{34,1}(:,2))
        hold on
%         end
    hold on 
    
    %Calculate dwell time per quadrant 
    for i=1:size(xi,2)-1
    tempXin=[0 xi(i) xi(i+1)];
    tempYin=[0 yi(i) yi(i+1)];
    [in,~]=inpolygon(params.transcoords{34,1}(:,1),params.transcoords{34,1}(:,2),tempXin,tempYin);
    dwell(1,i)=sum(in)/fr;
    end
    
    xMax=params.Xmax(34); xMin=params.Xmin(34); yMax=params.Ymax(34); yMin=params.Ymin(34);
    
    xedge=linspace(-101,101,params.dia(34)/binsize);
    yedge=linspace(-101,101,params.dia(34)/binsize);
    
    [map] = histcounts2(params.transcoords{34,1}(:,1),params.transcoords{34,1}(:,2),xedge,yedge);
    map=flipud(imrotate(map,90));
    
    % smooth with gaussian and sigma of 1.5
    map = imgaussfilt(map, 1.5);
    
    % remove nan
%     map(map==0)=NaN;
map=params.upHBmap{19,1};
    map(map==0)=NaN;
    figure
    imAlpha=ones(size(map));
    imAlpha(isnan(map))=0;
    imagesc(map,'AlphaData',imAlpha);
    axis xy; axis off; hold on; box off; axis image;
    hold on
    plot(rescale(x,1,size(map,2)),rescale(y,1,size(map,1)),'k')
    
    
    %%plots 2 smaller circle annuli
%     plot((sin(0:2*pi/1000:2*pi)*((params.dia(i)/2)-subtraction)),cos(0:2*pi/1000:2*pi)*((params.dia(i)/2)-subtraction))
%     plot((sin(0:2*pi/1000:2*pi)*((params.dia(i)/2)-(subtraction*2))),cos(0:2*pi/1000:2*pi)*((params.dia(i)/2)-(subtraction*2)))
%     hold on
%     plot((params.transcoords{i,1}(:,1)), params.transcoords{i,1}(:,2))
%     axis image
        