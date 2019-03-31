function [out]=bordermodulation(siz,ts,x,y,samplerate,spkbin,binsize)
% Border Modulation [Peyrache, Schieferstein, Buzsaki, 2017, DOI:10.1038/s41467-017-01908-3]
%
% Inputs:
%           siz: size of ratemap in bins
%           ts: timestamps with spike timestamps spliced in
%           x: x coordinates with spike position spliced in
%           y: y coordinates with spike position spliced in
%           samplerate: sample rate in hz
%           spkbin: binary indicating which ts included a spike
%           binsize: size of each bin in ratemap (cm) 
%
% Output:
%           out: regression coefs for [top, right, bottom, left]
%
% Ryan E Harvey
%
warning off
% innerboundary (15cm near the boundary)
innerwidth=15/binsize;

% rescale x and y to length of ratemap
x=rescale(x,1,siz);
y=rescale(y,1,siz);

% #1 CALC IFR
% bin over 50ms
nframes=round(samplerate*.05);
spikets=ts(spkbin);
videots=ts(~spkbin);
edges=videots(1:nframes:end);
[N,~] = histcounts(spikets,edges);
% convert to spikes/sec
N=N*(samplerate/nframes);
% smooth over 200ms
nframes=round(samplerate*.2);
ifr=smooth(N,nframes);

% #2 FIND XY FRAMES NEAR THE BORDER
x=x(ismember(unique(ts),edges));
y=y(ismember(unique(ts),edges));

top=[[innerwidth+1,siz];[siz-innerwidth,siz];[siz-innerwidth,siz-innerwidth];[innerwidth+1,siz-innerwidth];[innerwidth+1,siz]];
right=[[siz-innerwidth,siz-innerwidth];[siz,siz-innerwidth];[siz,innerwidth+1];[siz-innerwidth,innerwidth+1];[siz-innerwidth,siz-innerwidth]];
bottom=[[innerwidth+1,innerwidth+1];[siz-innerwidth,innerwidth+1];[siz-innerwidth,1];[innerwidth+1,1];[innerwidth+1,innerwidth+1]];
left=[[1,siz-innerwidth];[innerwidth+1,siz-innerwidth];[innerwidth+1,innerwidth+1];[1,innerwidth+1];[1,siz-innerwidth]];

topin = inpolygon(x,y,top(:,1),top(:,2));
rightin = inpolygon(x,y,right(:,1),right(:,2));
bottomin = inpolygon(x,y,bottom(:,1),bottom(:,2));
leftin = inpolygon(x,y,left(:,1),left(:,2));

% #3 REGRESS IFR AND BINARY OF WHEN THE ANIMAL IS NEXT TO EACH BORDER
b=glmfit(ifr,topin(1:end-1),'poisson','link','log');
top_cor_coef=b(2);

b=glmfit(ifr,rightin(1:end-1),'poisson','link','log');
right_cor_coef=b(2);

b=glmfit(ifr,bottomin(1:end-1),'poisson','link','log');
bottom_cor_coef=b(2);

b=glmfit(ifr,leftin(1:end-1),'poisson','link','log');
left_cor_coef=b(2);

out=[top_cor_coef,right_cor_coef,bottom_cor_coef,left_cor_coef];

warning on
end