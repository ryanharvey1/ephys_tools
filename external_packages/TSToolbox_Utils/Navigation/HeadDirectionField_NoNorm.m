function [h,B,mu,p,kappa] = HeadDirectionField(tsa,ang,GoodRanges,varargin)

% USAGE
%     [h,d,mu,p,kappa] = HeadDirectionField(tsa,ang,GoodRanges)
%     
%     Inputs:
%     tsa: a ts object (typically a cell!)
%     ang: a tsd object of head orientation
%     GoodRanges: an intervalSet object defining the time of valid pos tracking
%     nbBins: number of angle bins (optionnal; default=50)

nbBins = 25;
if ~isempty(varargin)
    if isnumeric(varargin{1})
        nbBins=varargin{1};
    else
        error('Nb of bins must be numeric')
    end
end

tsa = Restrict(tsa,GoodRanges);
B = 2*pi*(0:1:nbBins-1)/nbBins;

if length(tsa)~=0
    
    ang = Restrict(ang,GoodRanges);
    angt = Restrict(ang,tsa);

    B = 2*pi*(0:1:nbBins-1)/nbBins;
    h0 = hist(Data(ang),B);
    h = hist(Data(angt),B);
    
    dt = 1/median(diff(Range(ang,'s')));
    %h = dt*h./h0;
    h(isnan(h))=0;

    h = [h h(1)];
    B = [B B(1)];

    [mu, kappa, p] = CircularMean(Data(angt));
else
    h = zeros(nbBins+1,0);
    mu = NaN;
    kappa = NaN;
    p = NaN;
end