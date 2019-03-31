function [h,v,B,mu] = HeadDirectionField(tsa1,tsa2,ang,GoodRanges,varargin)
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

tau = 100;
Q = MakeQfromS(tsdArray({tsa1;tsa2}),tau*10);
Q = Restrict(Q,GoodRanges);
dQ = Data(Q);
mx = max(dQ(:));

    angt = Restrict(ang,Q);

B = 2*pi*(0:1:nbBins-1)/nbBins;
da = Data(angt);
da = abs(repmat(da,[1 length(B)-1])-repmat(B(1:end-1),[size(da,1) 1]));
da(da>pi) = 2*pi-da(da>pi);
[dummy,minIxA] = min(da');
    
    corr = [];
    for b=1:length(B)
        
        h2 = hist2d(dQ(minIxA==b,1),dQ(minIxA==b,2),[0:mx],[0:mx]);
        %h2 = h2./repmat(sum(h2),[mx+1 1]);
        c = corrcoef(dQ(minIxA==b,1),dQ(minIxA==b,2));
        corr = [corr;c(1,2)];
%         figure(1),clf
%         imagesc(tanh(h2))
%         pause
    end
    
    keyboard
    
    m = h2'*[0:mx]';
    v = h2'*([0:mx].^2)' - m.^2;
    
    
    keyboard
    
    dt = 1/median(diff(Range(ang,'s')));
    h = dt*h./h0;
    h(isnan(h))=0;

    h = [h h(1)];
    B = [B B(1)];

    [mu, kappa, p] = CircularMean(Data(angt));
