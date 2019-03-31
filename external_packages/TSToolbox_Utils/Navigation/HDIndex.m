function [HdI,pHdI] = HDIndex(S,ang,ep,varargin)


% USAGE
%     HdI = HDIndex(S{c},ang,ep)
%     
%     Inputs:
%     S: spike train
%     ang: a tsd object of head orientation
%     ep: an intervalSet object defining the time of valid pos tracking
%
%     options:
%     [h,d,mu,p,kappa] = HDIndex(S,ang,ep,nbBins,sdSmooth)
%     nbBins: number of angle bins (default=360)
%     sdSmooth: s.d. of tuning curve smoothing (in number of bins, default=6)

%Adrien Peyrache, 2016

%Default values
nbBins = 60;
nbRnd = 1000;

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

S = Restrict(S,ep);
B = 2*pi*(0:1:nbBins-1)'/nbBins;

HdI = zeros(length(S),1);
pHdI = zeros(length(S),1);

ang = Restrict(ang,ep);
h0 = hist(mod(Data(ang),2*pi),B);
dt = 1/median(diff(Range(ang,'s')));

for c=1:length(S)
    if length(S{c})>100

        angt = Restrict(ang,S{c});    
        h = hist(mod(Data(angt),2*pi),B);
        h = dt*h./h0; 
        h(isnan(h))=0;

        HdI(c)  = norm(h(:)'*(cos(B)+1i*sin(B)))/sum(h(:));
    end
end
resVecShuf  = ones(length(S),nbRnd);
angD        = Data(ang);
angT        = Range(ang);
N           = length(angD);

for r=1:nbRnd

    angShuf = tsd(angT,angD(randperm(N)));
    for c=1:length(S)
    if length(S{c})>100
        angt    = Restrict(angShuf,S{c});
        h       = hist(mod(Data(angt),2*pi),B);
        h       = dt*h./h0;
        resVecShuf(c,r)  = norm(h(:)'*(cos(B)+1i*sin(B)))/sum(h(:));
    end
    end
end
pHdI = sum(resVecShuf>repmat(HdI,[1 nbRnd]),2)/nbRnd;
