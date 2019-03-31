function [x,y,fieldarea]=findFields2D(map,thres)

rates=sort(reshape(map,[],1),'descend');
rates(isnan(rates))=[];

maxarea=1;
i=1;
while i+1<length(rates) && maxarea<10
    clear x y z fieldarea
    tempmap=map;
    tempmap(map<rates(i)*thres)=0;
    
    C=contour(tempmap,[min(tempmap(tempmap>0)),min(tempmap(tempmap>0))]);
    
    if isempty(C)
        [r,c]=find(~isnan(map));
        k=convhull(c,r);
        x{1}=c(k);
        y{1}=r(k);
        fieldarea=polyarea(x{1},y{1})*3;
        break
    end
    i=i+1;
    
    
    m(1)=1;
    n=1;
    try
        while n<length(C)
            n=n+1;
            m(n) = m(n-1)+C(2,m(n-1))+1;
        end
    catch
    end
    
    for nn = 1:n-2
        x{nn} = C(1,m(nn)+1:m(nn+1)-1);
        y{nn} = C(2,m(nn)+1:m(nn+1)-1);
        z(nn) = C(1,m(nn));
    end
    
    for ii=1:length(x)
        fieldarea(ii)=polyarea(x{ii},y{ii})*3;
    end
    x(fieldarea<6)=[];
    y(fieldarea<6)=[];
    z(fieldarea<6)=[];
    
    maxarea=max(fieldarea);
    
    
end
% imagesc(map);colormap jet
% hold on% Allows plotting atop the preexisting peaks plot.
% for n = 1:length(x)
%     plot(x{n},y{n},'k', 'LineWidth', 2)
% end

% for n = find(z==-2) % now loop through the z = -2 values.
%     plot(x{n},y{n},'r:','linewidth',2)
% end