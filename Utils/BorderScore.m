function [b,dm]=BorderScore(RateMap)
% BorderScore- calculates border score and bins to wall
% Input:
%           RateMap
%
% Output: 
%           b: border score [-1 to 1] 1 indicates field is close to border 
%           dm: bins to wall
%
% Ryan Harvey

RateMap(isinf(RateMap))=NaN;

if sum(sum(isnan(RateMap)))==(size(RateMap,1)*size(RateMap,2))
    b=NaN;
    dm=NaN;
    return
end
if nansum(nansum(RateMap))==0
    b=NaN;
    dm=NaN;
    return
end

% FIND PEAK BIN & FIELD
Peak_Rate=max(max(RateMap)); r_20perc = Peak_Rate*0.30;

% FIND EDGES OF MAP
% [y_max,x_max]=size(RateMap);x_min=1;y_min=1;

% x=median([x_max,x_min]); y=median([y_max,y_min]); rad=(max(x_min-x_max,y_min-y_max))/2;
% th = 0:pi/179.5:2*pi;
% xunit = rad * cos(th) + x; yunit = rad * sin(th) + y;

[X,Y]=find(~isnan(RateMap));
k = convhull(X,Y);
xunit=X(k);
yunit=Y(k);

[y_20,x_20] = find(RateMap > r_20perc);
norm_RateMap = (RateMap - min(min(RateMap))) / ( max(max(RateMap)) - min(min(RateMap)) );

for j=1:length(x_20) % calculates min dist from every active bin to the walls
    for jj=1:length(xunit)
        d2(jj) = sqrt((x_20(j)-xunit(jj))^2+(y_20(j)-yunit(jj))^2);
    end
    d3(j)=min(d2)*norm_RateMap(y_20(j),x_20(j));
end
dm=mean(d3);


% FIND EDGE BINS
% Top
% for i=1:x_max
%     for ii=1:y_max
%         if ~isnan(RateMap(ii,i))
%             rt(i)=ii;
%             break
%         end
%     end
% end
% % Bottom
% for i=x_max:-1:1
%     for ii=y_max:-1:1
%         if ~isnan(RateMap(ii,i))
%             rb(i)=ii;
%             break
%         end
%     end
% end
% % concat
% edgeBins=[[[rt]';[fliplr(rb)]'],[[1:x_max]';[x_max:-1:1]']];
% % fill in sides
% edgeBins=[edgeBins;[[edgeBins(end,1)-1:-1:edgeBins(1,1)+1]',[ones(1,length(edgeBins(end,1)-1:-1:edgeBins(1,1)+1))]']];
% 
% edgeBins=[edgeBins(1:max(edgeBins),:);[[edgeBins(max(edgeBins(:,2)),1)+1:edgeBins(max(edgeBins(:,2)+1),1)-1]',...
%     [repelem(max(edgeBins(:,2)),length([edgeBins(max(edgeBins(:,2)),1)+1:edgeBins(max(edgeBins(:,2)+1),1)-1]))'];edgeBins(max(edgeBins):end,:)]];
% % sort
% [~,Index]=sort(edgeBins(:,2));
% edgeBins=edgeBins(Index,:);
% field=[y_20,x_20];

% [X,Y]=find(~isnan(RateMap));
% k = convhull(X,Y);

% figure;scatter(X,Y);hold on
% plot(X(k),Y(k),'r-')
% plot(x_20,y_20)

% scatter(edgeBins(:,1),edgeBins(:,2),'k')


cm=sum(ismember([x_20,y_20],[xunit,yunit],'rows'));

b=(cm-dm)/(cm+dm);

end