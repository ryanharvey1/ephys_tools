% RunDirandfieldshape
% 
% REVIEWER 1
% Minor Comment 2: p. 5, line 169. Is there an explanation for the elongated place fields?
% For example, does the elongation match the most frequent direction of travel?
% REVIEWER 3
% Comment 2: Behavior patterns in titled mice should be compared to the controls.
% If tilted mice tend to run in specific stereotypic patterns (for example, turn less),
% that might explain the elongation of place fields in titled mice in the center of the arena.
% The authors should check that this is not the case. Do titled mice run closer to the walls?
% Behavior should be analyzed more thoroughly for this data set.

clear;clc;close all
addpath ('/Users/ryanharvey/GoogleDrive/MatlabDir/CircStat2012a','/Users/ryanharvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis')
load('/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/NewOCC_Map_workspace2.mat')
FigureLocation='/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/Figures_For_Paper';
%% COLLECT PASSES THROUGH FIELD
[ResultsC]=analyze(ratemapC,dataC,0);
[ResultsT]=analyze(ratemapT,dataT,0);
ResultsT.cm2wall(94,1)=[];
ResultsT.perx(94,1)=[];
ResultsT.Ecc(94,1)=[];

% Plot results
[ AllStats ] = ScatterBox([ResultsC.perx,ResultsC.cm2wall,ResultsC.Ecc],[ResultsT.perx,ResultsT.cm2wall,ResultsT.Ecc],{'Control','Tilted'},{'Crossing Freq','Cm to Wall','Eccentricity'},2)

med=median([ResultsC.cm2wall;ResultsT.cm2wall]);
[ AllStats ] = ScatterBox([ResultsC.perx(ResultsC.cm2wall>med,:)],[ResultsT.perx(ResultsT.cm2wall>med,:)],{'Control','Tilted'},{'Away from Wall X freq'},2)
[ AllStats ] = ScatterBox([ResultsC.perx(ResultsC.cm2wall<med,:)],[ResultsT.perx(ResultsT.cm2wall<med,:)],{'Control','Tilted'},{'Close to Wall X freq'},2)

[ AllStats ] = ScatterBox([ResultsC.Ecc(ResultsC.cm2wall<med,:)],[ResultsT.Ecc(ResultsT.cm2wall<med,:)],{'Control','Tilted'},{'Close to Wall Eccentricity'},2)
[ AllStats ] = ScatterBox([ResultsC.Ecc(ResultsC.cm2wall>med,:)],[ResultsT.Ecc(ResultsT.cm2wall>med,:)],{'Control','Tilted'},{'Away from Wall Eccentricity'},2)

med=median([ResultsC.Ecc;ResultsT.Ecc]);med=med+.1;
[ AllStats ] = ScatterBox([ResultsC.perx(ResultsC.Ecc<med,:)],[ResultsT.perx(ResultsT.Ecc<med,:)],{'Control','Tilted'},{'Low Ecc'},2)
[ AllStats ] = ScatterBox([ResultsC.perx(ResultsC.Ecc'>med,:)],[ResultsT.perx(ResultsT.Ecc>med,:)],{'Control','Tilted'},{'High Ecc'},2)

print(figure(1),'-bestfit','-dpdf', '-r600',[FigureLocation,filesep,'Control_CELL18_PassField_EX.pdf'])

%% ////////////////////////////////////////////////////////////////////////
function [results]=analyze(ratemap,data,plots)
mapC=[];crossings=[];
sessions=fieldnames(ratemap);
sessions=sessions([1:5:length(sessions)],1);
% results=zeros(length(sessions),4);
for i=1:length(sessions)
    if i==94
        continue
    end
    disp([num2str(i),'  OF ',num2str(length(sessions)),'  CELLS  '])
    % pull in matrix
    mapC=ratemap.(sessions{i});
    [b,dm]=BorderScore(mapC);results.cm2wall(i,1)=dm*2.44;
    % find field
    [field,~]=FindFF2D(mapC);
    % find field
    [row,col]=find(field);
    %find boundary
    k=boundary(row,col);
    bound=[col(k),row(k)];
    results.Ecc(i,1)=Eccentricity(col(k),row(k));
    
    bound(:,1)=smooth(bound(:,1));
    bound(:,2)=smooth(bound(:,2));
    
    % find max diameter of field & calc angle
    for j=1:length(bound);for jj=1:length(bound);d2(jj)=sqrt((bound(j,1)-bound(jj,1))^2+(bound(j,2)-bound(jj,2))^2);end;d3(j)=max(d2);end;dia=find(d3==max(d3));clear d3 d2
    
    for t=2:length(dia)
        x=[bound(dia(1),1);bound(dia(t),1)];y=[bound(dia(1),2);bound(dia(t),2)];
        P1=[x(1,1),y(1,1)];
        P2=[x(end,1),y(end,1)];
        m=(P2(2)-P1(2))/(P2(1)-P1(1));
        N=100;
        x=linspace(0,50,N);
        y=P1(2)+m*(x-P1(1));
        x=x';y=y';
        P=interX([x,y]',bound');P=unique([round(P)]','rows')';
        if size(P,2)>1;d(t-1,:)=sqrt((P(1,1)-P(1,2))^2+(P(2,1)-P(2,2))^2);else;d(t-1)=0;end
    end
    [~,I]=max(d);dia=dia(1,[1,I+1]);
    
    x1=linspace(bound(dia(1),1),bound(dia(2),1),11);
    y1=linspace(bound(dia(1),2),bound(dia(2),2),11);
    
    theta=deg2rad(90);
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    rotData = R*([x1',y1']' - [median(x1),median(y1)]') + [median(x1),median(y1)]';
    x=rotData(1,:)';
    y=rotData(2,:)';
    
    P1=[x(1,1),y(1,1)];
    P2=[x(end,1),y(end,1)];
    % Extend line
    m = (P2(2)-P1(2))/(P2(1)-P1(1)) ;  % slope of line
    N = 100 ;
    x = linspace(0,30,N) ;  % specify your x limits
    y = P1(2)+m*(x-P1(1)) ;  % y = y1+m*(x-x1) ;
    x=x';
    y=y';
    P=interX([x,y]',bound'); P=unique([round(P)]','rows')';
    d=sqrt((P(1,1)-P(1,2))^2+(P(2,1)-P(2,2))^2);
    
    l1=2;l2=2;jj=.5;
    k1=1;k2=1;
    done1=0;
    done2=0;
    while (done1+done2)~=2
        if l1>1 && done1==0
            if sum(round(x)-round(y))==0
                P1=interX([x-jj,y+jj]',bound');
            else
                P1=interX([x+jj,y+jj]',bound');
            end
            P1=unique([round(P1)]','rows')';
            l1=size(P1,2);
            if l1>1 && done1==0
                d1=sqrt((P1(1,1)-P1(1,2))^2+(P1(2,1)-P1(2,2))^2);
                if l1>1 && d1>=d*.50
                    %                     plot(x-jj,y+jj,'LineWidth', 5, 'color', 'r');
                else
                    done1=1;
                    l1idx(k1)=jj-.7; k1=k1+1;
                end
            else
                done1=1;
                l1idx(k1)=jj-.7; k1=k1+1;
            end
        end
        if l2>1 && done2==0
            if sum(round(x)-round(y))==0
                P2=interX([x+jj,y-jj]',bound');
            else
                P2=interX([x-jj,y-jj]',bound');
            end
            P2=unique([round(P2)]','rows')';
            l2=size(P2,2);
            if l2>1 && done2==0
                d2=sqrt((P2(1,1)-P2(1,2))^2+(P2(2,1)-P2(2,2))^2);
                if l2>1 && d2>=d*.50
                    %                     plot(x+jj,y-jj,'LineWidth', 5, 'color', 'r');
                else
                    done2=1;
                    l2idx(k2)=jj-.7; k2=k2+1;
                end
            else
                done2=1;
                l2idx(k2)=jj-.7; k2=k2+1;
            end
        end
        jj=jj+.5;
    end
    if l1idx(1)<1 || l2idx(1)<1
        testflag=1
    end
    l1=[x+l1idx(1),y+l1idx(1)];
    l2=[x-l2idx(1),y-l2idx(1)];
    
    frames=data.(sessions{i});
    frames(:,2:3)=[smooth(rescale(frames(:,2),1,25),12),smooth(rescale(frames(:,3),1,25),12)];
    
    if plots==1
        fig=figure; fig.Color=[1 1 1]; set(fig,'Position',[1 5 720 800]);
        rad=12.2;th=0:pi/179.5:2*pi;
        xunit=rad*cos(th)+13;yunit=rad*sin(th)+13;
        % plot circle
        plot(xunit,yunit,'LineWidth', 5, 'color',[0 0 0]+.2);hold on
        % plot path
        plot(frames(:,2), frames(:,3), 'LineWidth', 3, 'color',[0,0,0]+.5);
        %         spk=frames(frames(:,5)==1,:);scatter(spk(:,2), spk(:,3), 70, 'filled', 'r');
    end
    
    in=inpolygon(frames(:,2),frames(:,3),bound(:,1),bound(:,2));
    
    % find durations field at least 200ms and over
    dsig = diff([0 (abs(in') >= eps) 0]);
    startIndex = find(dsig > 0);
    endIndex = find(dsig < 0)-1;
    stringIndex = (endIndex-startIndex+1 >= 12);
    startIndex = startIndex(stringIndex);
    endIndex = endIndex(stringIndex);
    indices = zeros(1,max(endIndex)+1);
    indices(startIndex) = 1;
    indices(endIndex+1) = indices(endIndex+1)-1;
    indices = find(cumsum(indices));
    in=zeros(length(in),1);
    in(indices',1)=1;
    
    pass=find(diff([0 indices])>1==1);
    for j=1:length(pass)
        if j+1>length(pass); break;end
        if pass==1;tempidx=indices; else; tempidx=indices(pass(j):pass(j+1)-1);end
        tempXY=[frames(tempidx',2),frames(tempidx',3)];
        %         if plots==1;plot(l1(:,1),l1(:,2),'LineWidth', 5, 'color', 'g');plot(l2(:,1),l2(:,2),'LineWidth', 5, 'color', 'g');end
        P1=interX(tempXY',[l1(:,1),l1(:,2)]');
        P1=unique([round(P1)]','rows')';
        lsize1=size(P1,2);
        P2=interX(tempXY',[l2(:,1),l2(:,2)]');
        P2=unique([round(P2)]','rows')';
        lsize2=size(P2,2);
        if lsize1>0 && lsize2>0
            crossings(j,1)=1;
            if plots==1
                plot(frames(tempidx',2),frames(tempidx',3),'LineWidth', 5, 'color',[0 0 0]+.2)
                %                 scatter(P1(1),P1(2), 60, 'filled', 'k');
                %                 scatter(P2(1),P2(2), 60, 'filled', 'k');
            end
        else
            crossings(j,1)=0;
        end
    end
    results.perx(i,1)=sum(crossings)/length(crossings);
    
    % plot field
    if plots==1;plot(bound(:,1),bound(:,2),'LineWidth', 6, 'color', 'r');axis off;axis image;box off;end
    if plots==1; pause(.0001);close all;end
    

    % %     ///////////////////////////testing/////////////////////
    %     pass=find(diff([0 indices])>1==1);
    %     for j=1:length(pass)
    %         if j+1>length(pass); break;end
    %         tempidx=indices(pass(j):pass(j+1)-1);
    % %         plot(frames(tempidx',2),frames(tempidx',3))
    %         x=linspace(min([bound(dia(1),1),bound(dia(2),1)]),max([bound(dia(1),1),bound(dia(2),1)]),length(tempidx));
    %         y=linspace(min([bound(dia(1),2),bound(dia(2),2)]),max([bound(dia(1),2),bound(dia(2),2)]),length(tempidx));
    %         r(j,1) = corr2([frames(tempidx',2),frames(tempidx',3)],[x',y']);
    %         normDeg=mod(frames(tempidx',7)-repmat(diaangle(1),length(tempidx),1),360);
    %         absDegDif=min(360-normDeg,normDeg);
    %         absDegDif(absDegDif>90)=absDegDif(absDegDif>90)-90;
    %         degerror(j,1)=mean(absDegDif);
    %     end

end
end
function P = interX(L1,varargin)
%INTERX Intersection of curves
%   P = INTERX(L1,L2) returns the intersection points of two curves L1
%   and L2. The curves L1,L2 can be either closed or open and are described
%   by two-row-matrices, where each row contains its x- and y- coordinates.
%   The intersection of groups of curves (e.g. contour lines, multiply
%   connected regions etc) can also be computed by separating them with a
%   column of NaNs as for example
%
%         L  = [x11 x12 x13 ... NaN x21 x22 x23 ...;
%               y11 y12 y13 ... NaN y21 y22 y23 ...]
%
%   P has the same structure as L1 and L2, and its rows correspond to the
%   x- and y- coordinates of the intersection points of L1 and L2. If no
%   intersections are found, the returned P is empty.
%
%   P = INTERX(L1) returns the self-intersection points of L1. To keep
%   the code simple, the points at which the curve is tangent to itself are
%   not included. P = INTERX(L1,L1) returns all the points of the curve
%   together with any self-intersection points.
%
%   Example:
%       t = linspace(0,2*pi);
%       r1 = sin(4*t)+2;  x1 = r1.*cos(t); y1 = r1.*sin(t);
%       r2 = sin(8*t)+2;  x2 = r2.*cos(t); y2 = r2.*sin(t);
%       P = InterX([x1;y1],[x2;y2]);
%       plot(x1,y1,x2,y2,P(1,:),P(2,:),'ro')

%   Author : NS
%   Version: 3.0, 21 Sept. 2010

%   Two words about the algorithm: Most of the code is self-explanatory.
%   The only trick lies in the calculation of C1 and C2. To be brief, this
%   is essentially the two-dimensional analog of the condition that needs
%   to be satisfied by a function F(x) that has a zero in the interval
%   [a,b], namely
%           F(a)*F(b) <= 0
%   C1 and C2 exactly do this for each segment of curves 1 and 2
%   respectively. If this condition is satisfied simultaneously for two
%   segments then we know that they will cross at some point.
%   Each factor of the 'C' arrays is essentially a matrix containing
%   the numerators of the signed distances between points of one curve
%   and line segments of the other.

%...Argument checks and assignment of L2
%     error(nargchk(1,2,nargin));
if nargin == 1,
    L2 = L1;    hF = @lt;   %...Avoid the inclusion of common points
else
    L2 = varargin{1}; hF = @le;
end

%...Preliminary stuff
x1  = L1(1,:)';  x2 = L2(1,:);
y1  = L1(2,:)';  y2 = L2(2,:);
dx1 = diff(x1); dy1 = diff(y1);
dx2 = diff(x2); dy2 = diff(y2);

%...Determine 'signed distances'
S1 = dx1.*y1(1:end-1) - dy1.*x1(1:end-1);
S2 = dx2.*y2(1:end-1) - dy2.*x2(1:end-1);

C1 = feval(hF,D(bsxfun(@times,dx1,y2)-bsxfun(@times,dy1,x2),S1),0);
C2 = feval(hF,D((bsxfun(@times,y1,dx2)-bsxfun(@times,x1,dy2))',S2'),0)';

%...Obtain the segments where an intersection is expected
[i,j] = find(C1 & C2);
if isempty(i),P = zeros(2,0);return; end;

%...Transpose and prepare for output
i=i'; dx2=dx2'; dy2=dy2'; S2 = S2';
L = dy2(j).*dx1(i) - dy1(i).*dx2(j);
i = i(L~=0); j=j(L~=0); L=L(L~=0);  %...Avoid divisions by 0

%...Solve system of eqs to get the common points
P = unique([dx2(j).*S1(i) - dx1(i).*S2(j), ...
    dy2(j).*S1(i) - dy1(i).*S2(j)]./[L L],'rows')';

    function u = D(x,y)
        u = bsxfun(@minus,x(:,1:end-1),y).*bsxfun(@minus,x(:,2:end),y);
    end
end
