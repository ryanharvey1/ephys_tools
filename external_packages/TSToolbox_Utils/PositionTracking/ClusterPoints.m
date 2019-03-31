function In = ClusterPoints(Data, IfPlot)
%function In = ClusterPoints(Data, IfPlot)
% Data is two column matrix of coordinates, IfPlot =1 if you want Data to
% be ploted by this function (0, if they are plotter already in current
% axes). 
% starts clustering manual interface like in klusters
% - left click - new line
%   right click - one back
%   middle click - close the polygon.
% don't click too fast -it wrongly interprets double click I guess.
%  returns In - binary vector where ones corresopond to caught fish :))


if IfPlot
    plot(Data(:,1), Data(:,2),'.','markersize',5);
end
linepos=1;
pl=[];%matrix of polyline coords
h_lastline=[];
hold on

while 1
    [x y button]=PointInput(1);

    switch button
        case 1 % left button
            pl(linepos,:)=[x y];
            if linepos>1
                h_lastline(end+1) = line([pl(linepos-1,1),pl(linepos,1)],[pl(linepos-1,2),pl(linepos,2)]);
                set(h_lastline(end),'Color','k');
            end;
            linepos=linepos+1;
        case 2% middle button
            if linepos>2
                pl(linepos,:)=pl(1,:);
                h_lastline(end+1) = line([pl(linepos-1,1),pl(linepos,1)],[pl(linepos-1,2),pl(linepos,2)]);
                set(h_lastline(end),'Color','k');
                break;
            end
        case 3 %right button
            if linepos>2
                linepos=linepos-1;
                delete(h_lastline(end));
                h_lastline(end)=[];
                pl(end,:)=[];

            end;
    end
end
In = inpolygon(Data(:,1), Data(:,2), pl(:,1),pl(:,2));
if nargout<1

    plot(Data(In,1),Data(In,2),'r*');
end


hold off;
