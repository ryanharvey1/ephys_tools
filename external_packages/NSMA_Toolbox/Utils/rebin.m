function [rx,ry] = rebin(x,y,nbins)
% 
%  [rx,ry] = rebin(x,y,nbins)
%
% rebin a histogram with x-bins x and contents y into 
% a histogram rx,ry with nbins < length(y), where
% rx is left edge of new bin and
% ry is averaged content over a block of nav=floor(length(y)/nbins)
% adjecent bins
%
% PL Dec.99


nb = length(y);
if nbins >= nb
   error('nbins larger than length(y)');
end

nav = floor(nb/nbins);
nrest = nb-(nbins-1)*nav;
if nav < nrest 
   nav = nav+1;
   nrest =nb-(nbins-1)*nav;
end
rx = zeros(nbins,1);
ry = rx;
for i=1:nbins-1
   rx(i) = x((i-1)*nav+1);
   ry(i) = sum(y((i-1)*nav+1:i*nav))/nav;
end
if nrest > 0
   rx(nbins) = x((nbins-1)*nav+1);
   ry(nbins) = sum(y((nbins-1)*nav+1:end))/nrest;
end
