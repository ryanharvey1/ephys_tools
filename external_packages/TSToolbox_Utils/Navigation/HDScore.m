function hdScore = HDScore(h,peak,B)


% USAGE
%     HDScore = HDScore(S,peak,B)

%Adrien Peyrache, 2016

hdScore = zeros(size(h,2),1);

for c=1:size(h,2)

    totFR   = sum(h(:,c));
    [~,ix]  = min(abs(B-peak(c)));
    ix      = ix(1);
    ii      = 0;
    ixA     = ix;
    
    while sum(h(ixA,c))<totFR/2
        ii      = ii+1;
        ixA     = mod(ix-ii:ix+ii,length(B))+1;
    end
    
    dAlpha = abs(B(ixA(end)) - B(ixA(1)));
    if dAlpha>pi
        dAlpha = 2*pi - dAlpha;
    end
    hdScore(c) = 1-dAlpha/pi;
end
