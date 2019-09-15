

range_=1:size(params.subID,1);

for j=range_
    if ~isempty(params.cueCoords{j})
        k=convhull(params.cueCM{j}(:,1),params.cueCM{j}(:,2));
    end
    fig=figure;fig.Color=[1 1 1];
    plot(sin(0:pi/360:2*pi)*(202/2),cos(0:pi/360:2*pi)*(202/2),'k')
    hold on; plot((cos(linspace(-pi,pi,1000))+0)*(101*.8),(sin(linspace(-pi,pi,1000))+0)*(101*.8),'--k')
    
    hold on
    plot(params.noseCM{j}(:,1),params.noseCM{j}(:,2),'r')
    plot(params.headCM{j}(:,1),params.headCM{j}(:,2),'b')
    plot(params.backCM{j}(:,1),params.backCM{j}(:,2),'g')
    plot(params.buttCM{j}(:,1),params.buttCM{j}(:,2),'c')
    hold on
    legend('maze','80% annulus','nose','head','back','butt','location','best')
    if ~isempty(params.cueCoords{j})
        plot(params.cueCM{j}(k,1),params.cueCM{j}(k,2),'-r','LineWidth',2)
    end
    axis square
    
end


%plot 3d paths
range_=1:size(params.subID,1);

for j=range_
    if ~isempty(params.cueCoords{j})
        k=convhull(params.cueCM{j}(:,1),params.cueCM{j}(:,2));
    end
    fig=figure;fig.Color=[1 1 1];
    plot(sin(0:pi/360:2*pi)*(202/2),cos(0:pi/360:2*pi)*(202/2),'k')
    hold on; plot((cos(linspace(-pi,pi,1000))+0)*(101*.8),(sin(linspace(-pi,pi,1000))+0)*(101*.8),'--k')
    
    hold on
    plot3(params.backCM{j}(:,1),params.backCM{j}(:,2),params.ts{j}(:,1),'k')
    hold on
    legend('back''location','best')
    if ~isempty(params.cueCoords{j})
        plot(params.cueCM{j}(k,1),params.cueCM{j}(k,2),'-r','LineWidth',2)
    end
    axis square
    
end