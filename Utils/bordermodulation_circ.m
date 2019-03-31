function [out]=bordermodulation_circ(ratemap,ts,x,y,samplerate,spkbin,binsize)
% Border Modulation for cylinder
% Adapted from [Peyrache, Schieferstein, Buzsaki, 2017, DOI:10.1038/s41467-017-01908-3]
%
%
%
% Algorithm:
% 1. Find the center of firing field
%       (contour based on 20% of peak rate & at least 10 bins)
% 2. Finds the angular difference between the center of your ratemap and the
%       center of your field & rotates xy points by this difference.
% 3. Calculates IFR over 100ms bins
% 4. Locates xy points at least 15cm from the arena's border
% 5. Splits up xy points into 4 90 degree slices
% 6. Regresses IFR to the occupancy in each quarter to find a relationship
%       between the enviroments boundaries and firing rate.
%
% Inputs:
%           ratemap: ratemap
%           ts: timestamps(sec) with spike timestamps spliced in
%           x: x coordinates with spike position spliced in
%           y: y coordinates with spike position spliced in
%           samplerate: sample rate in hz
%           spkbin: binary indicating which ts included a spike
%           binsize: size of each bin in ratemap (cm)
%
% Output:
%           out: regression coefs for [top, right, bottom, left]
%
%
% Note: Low spiking cells will have erroneous values
%
% Ryan E Harvey
%
warning off
% FIND CENTER OF FIRING FIELD AND ROTATE XY SO THAT FIELD WILL BE LINED
% UP WITH ZERO
siz=length(ratemap);
[field,~]=FindFF2D(ratemap);
[r,c]=find(field==1);

    try
    k=convhull(c,r);
    catch
        out=NaN;
        return
    end
    fieldbound=[c(k),r(k)];
    % find center of field
    [fcx,fcy] = centroid(polyshape(fieldbound(:,1),fieldbound(:,2)));
    % find out how much the center of the field deviates from 0 degrees
    degoffset=rad2deg(atan2(fcy-siz/2,fcx-siz/2));
    % rescale x and y to length of ratemap
    x=rescale(x,1,siz);
    y=rescale(y,1,siz);
    % Rotate xy coordinates
    % Create rotation matrix
    theta=-degoffset;
    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    % Rotate your point(s)
    point = [x y]';
    rotpoint = R*point;
    x=rotpoint(1,:);
    y=rotpoint(2,:);
    % locate center of rotated xy
    xtemp=x(~isnan(x));
    ytemp=y(~isnan(y));
    k=convhull(xtemp,ytemp);
    fieldbound=[xtemp(k)',ytemp(k)'];
    [cx,cy] = centroid(polyshape(fieldbound(:,1),fieldbound(:,2)));
    % INSIDE BOUND
    r=(siz-binsize*2)/2;
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + cx;
    yunit = r * sin(th) + cy;
    
    % CALC IFR
    % bin over 100ms
    nframes=round(samplerate*.100);
    spikets=ts(spkbin);
    videots=ts(~spkbin);
    
    edges=videots(1:nframes:end);
    [N,~] = histcounts(spikets,edges);
    % convert to spikes/sec
    N=N*(samplerate/nframes);
    % smooth over 200ms
    nframes=round(samplerate*.2);
    ifr=smooth(N,nframes);
    
    % FIND XY FRAMES NEAR THE BORDER
    x=x(ismember(unique(ts),edges));
    y=y(ismember(unique(ts),edges));
    
    % INSIDE BOUND
    in=inpolygon(x,y,xunit,yunit);
    donutidx=find(~in);
    
    % BOUND FOR EACH QUARTER OF ARENA
    right=[cx,cy;max(x),min(y);max(x),max(y);cx,cy];
    top=[cx,cy;max(x),max(y);min(x),max(y);cx,cy];
    left=[cx,cy;min(x),max(y);min(x),min(y);cx,cy];
    bottom=[cx,cy;min(x),min(y);max(x),min(y);cx,cy];
    
    % INDEX OF XY POINTS IN ABOVE QUARTERS
    %  [top, right, bottom, left]
    shapes{1}=find(inpolygon(x,y,top(:,1),top(:,2)));
    shapes{2}=find(inpolygon(x,y,right(:,1),right(:,2)));
    shapes{3}=find(inpolygon(x,y,bottom(:,1),bottom(:,2)));
    shapes{4}=find(inpolygon(x,y,left(:,1),left(:,2)));
    
    % REGRESS IFR AND BINARY OF WHEN THE ANIMAL IS NEXT TO EACH BORDER
    for i=1:length(shapes)
        inShape=zeros(length(ifr),1);
        idx=intersect(donutidx,shapes{i});
        idx(idx>length(inShape))=[];
        inShape(idx)=1;
        inShape=logical(inShape);
        b=glmfit(ifr,inShape,'poisson','link','log');
        out(1,i)=b(2);
    end

warning on
end