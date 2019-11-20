function [BW,maskedImage,x,y,fieldarea,X] = segmentImage(varargin)
%segmentImage Segment image using kmeans clustering 
%  [BW,maskedImage,x,y,fieldarea,X] = segmentImage(varargin) segments image X
%  
%  
%  Input:  
%           map: map to segment (default peaks)
%           binsize: size of each bin (default 3)
%           figs: do you want figures? (default false)
%           upscalefac: up scale factor (default 15). Increases image size which
%               effectively smooths image which helps segmentation. 
%
%  Output:
%           BW: logical of segmentation
%           maskedimage: map with segmentaion applied 
%           x: x values for segmentation boundaries (cell array)
%           Y: y values for segmentation boundaries (cell array)
%           fieldarea: segmentation area
%           X: map with nans accounted for 
% 
% Ryan Harvey
%----------------------------------------------------

% parse input
p = inputParser;
p.addParameter('map',peaks);
p.addParameter('binsize',3);
p.addParameter('figs',false);
p.addParameter('upscalefac',15);
p.parse(varargin{:});

X = p.Results.map;
binsize = p.Results.binsize;
figs = p.Results.figs;
upscalefac = p.Results.upscalefac;

warning off

nanloc=padarray(isnan(X),[1 1],1,'both');
nanloc=imresize(nanloc,size(nanloc)*upscalefac);

X=padarray(X,[1 1],0,'both');

X(isnan(X))=0;
X=imresize(X,size(X)*upscalefac);

% Normalize input data to range in [0,1].
Xmin = min(X(:));
Xmax = max(X(:));
if isequal(Xmax,Xmin)
    X = 0*X;
else
    X = (X - Xmin) ./ (Xmax - Xmin);
end

% Auto clustering
sz = size(X);
im = single(reshape(X,sz(1)*sz(2),[]));
im = im - mean(im);
im = im ./ std(im);
s = rng;
rng('default');
L = kmeans(im,2,'Replicates',5);
rng(s);
BW = L == 2;
BW = reshape(BW,[sz(1) sz(2)]);

% Create masked image.
maskedImage = X;
maskedImage(~BW) = 0;

[C]=contourc(double(~BW),[1 1]);

if isempty(C)
    [r,c]=find(~isnan(~BW));
    k=convhull(c,r);
    x{1}=c(k);
    y{1}=r(k);
    fieldarea=polyarea(x{1},y{1})*binsize;
    if figs
        imagesc(X);hold on
        axis off
        axis image
        colormap jet
        plot(x{1},y{1},'r','LineWidth', 2)
    end
    return
end

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
    fieldarea(nn)=(polyarea(x{nn},y{nn})/upscalefac^2)*binsize;
end
% for some reason contourc gives multiple similar countours, so here we
% find and remove the duplicates within a certain tol.
% This is not the best way to do this
% [~,I]=uniquetol(cellfun('length',x),1);
%  zz=zeros(length(fieldarea),1);
% zz(I)=1;

% 
% x(fieldarea<30 | ~zz')=[];
% y(fieldarea<30 | ~zz')=[];
% fieldarea(fieldarea<30 | ~zz')=[];

X(nanloc)=NaN;
if figs
    imAlpha=ones(size(X));
    imAlpha(isnan(X))=0;
    imagesc(X,'AlphaData',imAlpha);
    axis off
    axis image
    colormap(viridis(255))
    hold on
    for n = 1:length(x)
        plot(x{n},y{n},'r','LineWidth', 3)
    end
end
warning off
end