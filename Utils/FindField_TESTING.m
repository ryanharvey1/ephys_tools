% FindField_TESTING
clear 
close all
figure

load('/Users/ryanharvey/Downloads/ratemaps.mat')
ratemap=ratemap(:,3:4);



[p,n]=numSubplots(20);
for per=1:20
% map=peaks(30);
% fornan=rand(30)>.95;
maptrue=ratemap{per,2};
map=ratemap{per,2};
map(isnan(map))=0;
% filtWidth = [2,2]; filtSigma = 1;
% imageFilter=fspecial('gaussian',filtWidth,filtSigma);
% ranconv=conv2(rand(size(map)),imageFilter,'same');
% maptrue=map.*ranconv;

% map(map<max(map(:))*.2)=.001;
% filtWidth = [3,3]; filtSigma = 2;
% imageFilter=fspecial('gaussian',filtWidth*filtSigma,filtSigma);
% imfilter(map,)

% SMOOTH
h = 1;
myfilter = fspecial('gaussian',[4 4]*h, h);
map = imfilter(map,myfilter,'replicate');
% map=conv2(maptrue,imageFilter,'same');
map(isnan(maptrue))=NaN;

% maptrue(fornan)=NaN;
% 

%  [P,logi]=FastPeakFind(map,max(map(:))*.2,imageFilter,2);
%  X=P(1:2:end);
%  Y=P(2:2:end);
%  peakval=map(logical(logi));
   
 
subplot(p(1),p(2),per)
[x,y,fieldarea]=findFields2D(map,.2);
pcolor(maptrue);hold on
for n = 1:length(x)
    plot(x{n},y{n},'k', 'LineWidth', 2)
end
% plot(X,Y,'*k')
axis image
colormap jet
shading flat
box off
axis off

fields{per}.x=x;
fields{per}.y=y;
fields{per}.fieldarea=fieldarea;


  
% subplot(4,3,per)
% [B,L] = bwboundaries(map>max(map(:))*.2,'noholes');
% pcolor(maptrue)
% hold on
% for k = 1:length(B)
%    bound = B{k};
%    plot(bound(:,2), bound(:,1), 'k', 'LineWidth', 2)
% end
% axis image
% colormap jet
% shading flat
% box off
% axis off
% title('bwboundaries method')
% 
% % subplot(4,3,per+1)
% [x,y,fieldarea]=findFields2D(inpaint_nans(map));
% pcolor(maptrue);hold on
% for n = 1:length(x)
%     plot(x{n},y{n},'k', 'LineWidth', 2)
% end
% axis image
% colormap jet
% shading flat
% box off
% axis off
% title('contours method')

% subplot(4,3,per+2)
% [field,fieldWidth]=FindFF2D(map);
% pcolor(maptrue);hold on
% [row,col]=find(field);
% k=boundary(row,col);
% plot(col(k),row(k),'k', 'LineWidth', 2)
% axis image
% colormap jet
% shading flat
% box off
% axis off
% title('Ryan Laura method')
end
