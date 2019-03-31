function Phasemap(ts,x,spkts,lfp,bins,trodeID,samplerate)

xedge=linspace(min(x),max(x),bins+1);
phaseedge=linspace(0,720,bins+1);

xspk=interp1(ts,x,spkts);

phase=interp1(lfp.ts,...
    lfp.theta_phase(trodeID,:),spkts','linear');
spkmap=histcounts2([xspk;xspk],[phase';phase'+2*pi]*180/pi,xedge,phaseedge);

phase=interp1(lfp.ts,...
    lfp.theta_phase(trodeID,:),ts','linear');
occ=histcounts2([x;x],[phase';phase'+2*pi]*180/pi,xedge,phaseedge);
occ=occ/samplerate;


phasemap=spkmap./occ;

phasemap(isnan(phasemap)) = 0;
phasemap(isinf(phasemap)) = 0;

% SMOOTH
h = 1.5;
h=4;
myfilter = fspecial('gaussian',[4 12]*h, h);
phasemap = imfilter([phasemap,phasemap,phasemap],myfilter,'replicate');
phasemap=phasemap(:,bins:(bins-1)*2);

pcolor(flipud(rot90(phasemap)));shading flat;box off;axis off;axis tight;
end
% filtWidth = [7 7]; filtSigma = 7;
% imageFilter=fspecial('gaussian',filtWidth,filtSigma);
% phasemap=conv2([phasemap,phasemap,phasemap],imageFilter,'same');
% phasemap=phasemap(:,bins:(bins-1)*2);
% pcolor(flipud(rot90(phasemap)));shading flat;box off;axis off;axis tight;