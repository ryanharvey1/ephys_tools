function out=egocentricmodulation(x,y,a,framesEXP,siz,samplerate,binsize)
% Egocentric Modulation for square enviorment [Peyrache, Schieferstein, Buzsaki, 2017, DOI:10.1038/s41467-017-01908-3]
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

% find innerboundary (15cm near the boundary)
innerwidth=15/binsize;
innerboundary=[[innerwidth+1,siz-innerwidth];[siz-innerwidth,siz-innerwidth];...
    [siz-innerwidth,innerwidth+1];[innerwidth+1,innerwidth+1];[innerwidth+1,siz-innerwidth]];
% figure;plot(x,y);hold on
% plot(innerboundary(:,1),innerboundary(:,2))

% in = inpolygon(x,y,innerboundary(:,1),innerboundary(:,2));
% xinzone=x(~in);
% yinzone=y(~in);

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

% boundary=[[linspace(1,siz,100)',repmat(siz,1,100)'];...
%     [repmat(siz,1,100)',linspace(siz,1,100)'];...
%     [linspace(siz,1,100)',ones(1,100)'];...
%     [ones(1,100)',linspace(1,siz,100)']];

boundary=[[(1:siz)',repmat(siz,1,siz)'];...
    [repmat(siz,1,siz)',linspace(siz,1,siz)'];...
    [linspace(siz,1,siz)',ones(1,siz)'];...
    [ones(1,siz)',(1:siz)']];


% top=[linspace(1,siz,100)',repmat(siz,1,100)'];
% right_=[repmat(siz,1,100)',linspace(siz,1,100)'];
% bottom=[linspace(siz,1,100)',ones(1,100)'];
% left_=[ones(1,100)',linspace(1,siz,100)'];

top=[linspace(1,siz,siz)',repmat(siz,1,siz)'];
right_=[repmat(siz,1,siz)',linspace(siz,1,siz)'];
bottom=[linspace(siz,1,siz)',ones(1,siz)'];
left_=[ones(1,siz)',linspace(1,siz,siz)'];

rightout=zeros(length(xinzone),1);
leftout=zeros(length(yinzone),1);
% d2=zeros(1,length(boundary));

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
    
    distances = sqrt(sum(bsxfun(@minus, boundary, [xinzone(j),yinzone(j)]).^2,2));
    closestwall = boundary(distances==min(distances),:);
    
    % check which wall is the closest to the rat  
    n=find([ismember(closestwall,top,'rows'),ismember(closestwall,right_,'rows'),...
        ismember(closestwall,bottom,'rows'),ismember(closestwall,left_,'rows')]);
        
    if length(n)>1
        n=n(randperm(length(n)));
    end
    switch n(1)
        case 1
            wall=top;
        case 2
            wall=right_;
        case 3
            wall=bottom;
        case 4
            wall=left_;
    end
    
    angle1 = angles(j)+60;
    angle2 = angles(j)-60;
    
    x2 = xinzone(j) + siz * cosd(angle1);
    y2 = yinzone(j) + siz * sind(angle1);
    
    x3 = xinzone(j) + siz * cosd(angle2);
    y3 = yinzone(j) + siz * sind(angle2);
    
    Pr = isintersect([[xinzone(j),yinzone(j)];[x2,y2]], [wall(1,:);wall(end,:)]);
    Pl = isintersect([[xinzone(j),yinzone(j)];[x3,y3]], [wall(1,:);wall(end,:)]);
    
    if plots==1
        if Pr
            [x_intersect,y_intersect]=intersectpoint([[xinzone(j),yinzone(j)];[x2,y2]],[wall(1,:);wall(end,:)]);
            scatter(x_intersect,y_intersect,'*y')
        end
        if Pl
            [x_intersect,y_intersect]=intersectpoint([[xinzone(j),yinzone(j)];[x3,y3]],[wall(1,:);wall(end,:)]);
            scatter(x_intersect,y_intersect,'*r')
        end
        p.XData = [xinzone(j);x2];
        p.YData = [yinzone(j);y2];
        p2.XData = xinzone(j);
        p2.YData = yinzone(j);
        p3.XData = [xinzone(j);x3];
        p3.YData = [yinzone(j);y3];
        drawnow
    end
    
    rightout(j,1)=Pl;
    leftout(j,1)=Pr;
    
    if Pr == Pl
        leftout(j,1)=0;
        rightout(j,1)=0;
    end
end
b=glmfit(ifr,leftout(1:end-1),'poisson','link','log');
left_cor_coef=b(2);

b=glmfit(ifr,rightout(1:end-1),'poisson','link','log');
right_cor_coef=b(2);

out=[left_cor_coef,right_cor_coef];

warning on
end

function intersect = isintersect(line1, line2)
A = [line1(1,:) - line1(2,:); line2(2,:) - line2(1,:)]';
if det(A) == 0
    intersect = 0;
else
    mu = A \ (line2(2,:) - line1(2,:))';
    intersect = all(mu >= 0) && all(mu <= 1);
end
end

function [x_intersect,y_intersect]=intersectpoint(line1, line2)
p1 = polyfit(line1(:,1),line1(:,2),1);
p2 = polyfit(line2(:,1),line2(:,2),1);
x_intersect = fzero(@(x) polyval(p1-p2,x),3);
y_intersect = polyval(p1,x_intersect);
end
