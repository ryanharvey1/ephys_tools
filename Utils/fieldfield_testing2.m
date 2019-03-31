

clear 
close all
figure

load('/Users/ryanharvey/Downloads/ratemaps.mat')
ratemap=ratemap(:,3:4);
[p,n]=numSubplots(20);
for per=1:20
maptrue=ratemap{per,2};

h = 1.5;
myfilter = fspecial('gaussian',[4 4]*h, h);
map = imfilter(maptrue,myfilter,'replicate');
% map=conv2(maptrue,imageFilter,'same');
map(isnan(maptrue))=NaN;
map=padarray(maptrue,[1 1],0,'both');

[zmax,imax,zmin,imin] = extrema2(map);
[I,J] = ind2sub(size(map),imax);
% figure;
% surf(map);hold on;colormap jet
% plot3(J,I,zmax,'*k')
% figure;
[~,logi]=FastPeakFind(map);
% imagesc(map);hold on
% plot(p(1:2:end),p(2:2:end),'r+')
 zmax=map(logical(logi));

map(isnan(map))=0;

for i=1:length(zmax)

    C=contour(map,[zmax(i)*.2,zmax(i)]);
    
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
        z{nn} = C(1,m(nn));
    end
%     figure
    for nn=1:length(x)
%         plot(x{nn},y{nn});hold on
%         plot(J(i),I(i),'*r');hold on
        infield(nn)=inpolygon(J(i),I(i),x{nn},y{nn});
    end
    
    x=x{infield};
    y=y{infield};
    z=z{infield};

    fieldarea(i)=polyarea(x,y)*3;
        
    X{i}=x;
    Y{i}=y;     
    clear x y z
end
subplot(p(1),p(2),per)
pcolor(maptrue);hold on
for n = 1:length(X)
    plot(X{n},Y{n},'k', 'LineWidth', 2)
end
end


