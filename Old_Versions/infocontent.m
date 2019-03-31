function out=infocontent(SmoothRateMap,occ)
[nBinsx,nBinsy]=size(SmoothRateMap);
rY=reshape(SmoothRateMap,nBinsx*nBinsy,1);
rY(isnan(rY) | isinf(rY))=0;
occRSHP=reshape(occ,nBinsx*nBinsy,1);
occSUM=sum(occRSHP);
pX=occRSHP./occSUM;
[nBins,nCells]=size(rY);
relR=rY./kron(ones(nBins,1),pX'*rY);
log_relR=log2(relR);
log_relR(isinf(log_relR))=0;
out=sum(kron(pX,ones(1,nCells)).*relR.*log_relR);
end