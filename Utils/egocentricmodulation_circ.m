function out=egocentricmodulation_circ(x,y,a,framesEXP,siz,samplerate,binsize)
% Egocentric Modulation for circular enviorment 
% based on [Peyrache, Schieferstein, Buzsaki, 2017, DOI:10.1038/s41467-017-01908-3]
% Inputs:
%           x: x coordinates with spike position spliced in
%           y: y coordinates with spike position spliced in
%           a: angles with spike angles spliced in
%           framesEXP: [ts,x,y,a,[],spk binary]
%           siz: size of ratemap in bins
%           samplerate: sample rate in hz
%           binsize: size of each bin in ratemap (cm) 
%
% Output:
%           out: regression coefs for [left side of body, right side of body]
%
% Ryan E Harvey
%
warning off

plots=0;

if ~any(a>pi | a>2*pi)
    a=rad2deg(a);
end

% rescale x and y to length of ratemap
x=rescale(x,1,siz);
y=rescale(y,1,siz);

% find 15cm from wall
r=(siz-15/binsize)/2;
th = 0:pi/50:2*pi;
xunit = r * cos(th) + median(1:siz);
yunit = r * sin(th) + median(1:siz);

innerboundary=[xunit', yunit'];

% bin over 50ms
nframes=round(samplerate*.05);
spikets=framesEXP(framesEXP(:,6)==1, 1);
videots=framesEXP(framesEXP(:,6)==0,1);
edges=videots(1:nframes:end);
[N,~] = histcounts(spikets,edges);
% convert to spikes/sec
N=N*(samplerate/nframes);
% smooth over 200ms
nframes=round(samplerate*.2);
ifr=smooth(N,nframes);

xinzone=x(ismember(unique(framesEXP(:,1)),edges));
yinzone=y(ismember(unique(framesEXP(:,1)),edges));
angles=a(ismember(unique(framesEXP(:,1)),edges));

r=siz/2;
th = 0:pi/50:2*pi;
xunit = r * cos(th) + median(1:siz);
yunit = r * sin(th) + median(1:siz);

boundary=[xunit', yunit'];

rightout=zeros(length(xinzone),1);
leftout=zeros(length(yinzone),1);

if plots==1
    close all
    fig=figure;fig.Color=[1 1 1];
    plot(boundary(:,1),boundary(:,2),'k');hold on
    box off;axis off;axis image
    p=plot(0,0,'y');
    p2=plot(0,0,'*k');
    p3=plot(0,0,'r');
    xlim([1-1 siz+1])
    ylim([1-1 siz+1])
end

for j=1:length(xinzone) % calculates min dist from every active bin to the walls
    if inpolygon(xinzone(j),yinzone(j),innerboundary(:,1),innerboundary(:,2)) || isnan(xinzone(j))
        continue
    end
    
    angle1 = angles(j)+60;
    angle2 = angles(j)-60;
    
    x2 = xinzone(j) + siz * cosd(angle1);
    y2 = yinzone(j) + siz * sind(angle1);
    
    x3 = xinzone(j) + siz * cosd(angle2);
    y3 = yinzone(j) + siz * sind(angle2);

    [I1,I2]=linecirc_intersect([xinzone(j),yinzone(j)],[x2,y2],siz/2,siz/2);
    [I3,I4]=linecirc_intersect([xinzone(j),yinzone(j)],[x3,y3],siz/2,siz/2);
    intersectz=[I1;I2;I3;I4];    
    distances = sqrt(sum(bsxfun(@minus,intersectz,[xinzone(j),yinzone(j)]).^2,2));
    [~,I]=min(distances);
    
    leftout(j,1)=I<3;
    rightout(j,1)=I>2;
    
    if plots==1
        if I<3
            scatter(intersectz(I,1),intersectz(I,2),'*b')
        else
            scatter(intersectz(I,1),intersectz(I,2),'*r')
        end
        p.XData = [xinzone(j);x2];
        p.YData = [yinzone(j);y2];
        p2.XData = xinzone(j);
        p2.YData = yinzone(j);
        p3.XData = [xinzone(j);x3];
        p3.YData = [yinzone(j);y3];
        drawnow
    end
end
b=glmfit(ifr,leftout(1:end-1),'poisson','link','log');
left_cor_coef=b(2);

b=glmfit(ifr,rightout(1:end-1),'poisson','link','log');
right_cor_coef=b(2);

out=[left_cor_coef,right_cor_coef];

warning on
end

function [I1,I2]=linecirc_intersect(P1,P2,C,r)
   A = P1-C;
   B = P2-P1;
   d2 = dot(B,B);
   t = (r^2-dot(A,A))*d2+dot(A,B)^2;
   Q = P1-dot(A,B)/d2*B;
   t2 = sqrt(t)/d2*B;
   I1 = Q + t2;
   I2 = Q - t2;
end
