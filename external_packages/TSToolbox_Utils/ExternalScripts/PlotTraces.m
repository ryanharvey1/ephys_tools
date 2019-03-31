%function h = PlotTraces(x,t,sr,scale,color,ifnorm, Axes2Fit,plotfun)
% x is a matrix of signals to plot: time x nchannels x ntraces 
% t is time axis, sr - sampling, scale - coefficient to scale up the signals
% color - to use, norm - if you want ot normolize them all evenly
function h = PlotTraces(x,varargin)
x=squeeze(x);
if ndims(x)<=2
    nChannels = min(size(x));
    nTime = max(size(x));
    if (size(x,1)~=nTime)
        x=x';
    end
    nTraces = 1;
else
    [nTime nChannels nTraces] = size(x);
end

[ t,sr,scale,color,ifnorm,Axes2Fit,plotfun] = DefaultArgs(varargin, { [], 1250,  1, 'k', 0, ylim,'plot'});
if isempty(t)
    t= [1:nTime]*1000/sr;
end
ChAmp = Srange(x);

%x = x*scale;
if nTraces ==1
    newx = x - repmat(mean(x([1:2,end-1:end],:),1),nTime,1);
    multipl = max(max(ChAmp),1);
    shift = diff(Axes2Fit)/(nChannels);
    newx = newx*scale./multipl*shift + shift + repmat([0:(nChannels-1)]*shift,nTime,1);
else
    
    newx = x - repmat(mean(x([1:2,end-1:end],:,:)),[nTime,1,1]);
    shift = max(median(sq(ChAmp),2));
    newx = newx - shift/2-repmat([0:nChannels-1]*shift,[nTime,1,nTraces]);
    newx = reshape(newx,nTime,[]);
end
%maxx = max(newx(:)); minx =min(newx(:));
%newx = Axes2Fit(2) - (newx-minx)*diff(Axes2Fit)/(maxx-minx);

if strcmp(plotfun,'plot')
    h= plot(t,newx,color);
else
    h= feval(plotfun,t,newx);
end
axis tight