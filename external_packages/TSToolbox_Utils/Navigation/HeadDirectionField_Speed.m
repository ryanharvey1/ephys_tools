function [h2,Bt,Bs] = HeadDirectionField_Speed(tsa1,ang,GoodRanges,varargin)
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

tau = 25.6;
Q = MakeQfromS(tsdArray({tsa1}),tau*10);
Q = Restrict(Q,GoodRanges);
dQ = Data(Q);
mx = max(dQ(:));

B = 2*pi*(0:1:nbBins-1)/nbBins;
    
    da = diff(Data(ang));
    da(da<pi) = 2*pi+da(da<pi);
    da(da>pi) = 2*pi - da(da>pi);
    da = gausssmooth(da,3,1);
    da = tsd(Range(ang),abs([da;0]));
    angSpd = Restrict(da,Q);
    
    angt = Restrict(ang,Q);
    da = Data(angt);
        
    Bt = 2*pi*(0:1:nbBins-1)/nbBins;
    mxSp = max(Data(angSpd));
    Bs = [0:0.01:0.25];
    
    h2 = histCont2d(dQ,Data(angt),Data(angSpd),Bt,Bs)';
    h2 = [h2;h2;h2];
    %h2 = gausssmooth(h2,1,1);
    h2 = h2(length(Bt)+1:2*length(Bt),:);
    h2 = h2/(tau/1000);
    
    