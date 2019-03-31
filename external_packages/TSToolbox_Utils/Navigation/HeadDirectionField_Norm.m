function [h,B,mu,p,kappa,h0] = HeadDirectionField(tsa,ang,GoodRanges,varargin)

% USAGE
%     [h,B,mu,p,kappa] = HeadDirectionField(tsa,ang,GoodRanges)
%     
%     Inputs:
%     tsa: a ts object (typically a cell!)
%     ang: a tsd object of head orientation
%     GoodRanges: an intervalSet object defining the time of valid pos tracking
%
%     options:
%     [h,d,mu,p,kappa] = HeadDirectionField(tsa,ang,GoodRanges,nbBins,sdSmooth)
%     nbBins: number of angle bins (default=360)
%     sdSmooth: s.d. of tuning curve smoothing (in number of bins, default=6)



%Default values
nbBins = 60;
sdSmooth = 1; %s.d. of smoothing filter in number of bins

h0 = [];
if ~isempty(varargin)
    if isnumeric(varargin{1}) && length(varargin{1})==1
        nbBins=varargin{1};
    else
        error('Nb of bins must be numeric')
    end
    if length(varargin) == 2
        if isnumeric(varargin{2}) && length(varargin{2})==1
            sdSmooth=varargin{2};
        else	
            error('S.d. of smoothing filter must be numeric')
        end
    end
end

tsa = Restrict(tsa,GoodRanges);
B = 2*pi*(0:1:nbBins-1)'/nbBins;
dt = 1/median(diff(Range(ang,'s')));

if ~isempty(Range(tsa))
    
    ang = Restrict(ang,GoodRanges);
    angt = Restrict(ang,tsa);
   
    h0 = hist(mod(Data(ang),2*pi),B);
    h = hist(mod(Data(angt),2*pi),B);
    
    h = dt*h./h0;
    
    h(isnan(h))=0;
    h = h(:);

    if sdSmooth
        h = gaussFiltAng(h,sdSmooth,1);
    end
    
    [mu, ~, p] = CircularMean(Data(angt));
    
    if mu>2*pi;mu = mu - 2*pi;end
    b = B;
    b(b>2*pi) = b(b>2*pi) - 2*pi;
    beta    = nlinfit(mu+B,h/sum(h)/median(diff(B)),'von_mises_pdf_params',[0 1]);
    kappa   = abs(beta(2));
    mu      = beta(1);
    h = [h;h(1)];
    B = [B;B(1)];

    [mu, ~, p] = CircularMean(Data(angt));

else
    h = zeros(nbBins+1,1);
    mu = NaN;
    kappa = NaN;
    p = NaN;
end