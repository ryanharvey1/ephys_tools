

range_=1:size(params.subID,1);

for i=range_
    if ~isempty(params.cueCoords{i})
        k=convhull(params.cueCM{i}(:,1),params.cueCM{i}(:,2));
    end
    fig=figure;fig.Color=[1 1 1];
    plot(sin(0:pi/360:2*pi)*(202/2),cos(0:pi/360:2*pi)*(202/2),'k')
    hold on; plot((cos(linspace(-pi,pi,1000))+0)*(101*.8),(sin(linspace(-pi,pi,1000))+0)*(101*.8),'--k')
    
    hold on
    plot(params.noseCM{i}(:,1),params.noseCM{i}(:,2),'r')
    plot(params.headCM{i}(:,1),params.headCM{i}(:,2),'b')
    plot(params.backCM{i}(:,1),params.backCM{i}(:,2),'g')
    plot(params.buttCM{i}(:,1),params.buttCM{i}(:,2),'c')
    hold on
    legend('maze','80% annulus','nose','head','back','butt','location','best')
    if ~isempty(params.cueCoords{i})
        plot(params.cueCM{i}(k,1),params.cueCM{i}(k,2),'-r','LineWidth',2)
    end
    axis square
    
end


