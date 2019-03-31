function plotheaddirectionmovie(x,y,theta)

siz=60;
fig=figure;fig.Color=[1 1 1];
plot(x, y,'.k');hold on
box off;axis off;axis image
p=plot(0,0,'r');
p2=plot(0,0,'*k');
p3=plot(0,0,'r');
xlim([min(x) max(x)])
ylim([min(y) max(y)])
for j=1:length(x)
    angle1 = theta(j)+60;
    angle2 = theta(j)-60;
    
    x2 = x(j) + siz * cosd(angle1);
    y2 = y(j) + siz * sind(angle1);
    
    x3 = x(j) + siz * cosd(angle2);
    y3 = y(j) + siz * sind(angle2);
    
    p.XData = [x(j);x2];
    p.YData = [y(j);y2];
    p2.XData = x(j);
    p2.YData = y(j);
    p3.XData = [x(j);x3];
    p3.YData = [y(j);y3];
    drawnow
end