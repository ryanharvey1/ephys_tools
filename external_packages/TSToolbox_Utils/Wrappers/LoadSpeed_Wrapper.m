function [linSpd,angSpd, GoodRanges] = LoadSpeed_Wrapper(fbasename,varargin)

% LoadSpeed_Wrapper loads linead and angular speed from whl or pos files
% USAGE:
% [linSpd,angSpd] = LoadSpeed_Wrapper(fbasename)

% Adrien Peyrache 2015


if ~isempty(varargin)
    w = varargin{1};
    if ~isa(w,'struct')
        error('argument must be a struct')
    end
    whl = w.whl;
    t = w.t;
    GoodRanges = w.GR;
else
    [whl,t,GoodRanges] = LoadPosition(fbasename);
end

if size(whl,2)~=4
    error('Need 2 LEDs !')
end

[folder dumy dumy]  = fileparts(fbasename);
if ~isempty(folder)
    folder =  [folder filesep];
end

if exist([folder 'Analysis/GeneralInfo.mat'],'file')
    warning off
    load([folder 'Analysis/GeneralInfo.mat'],'angOffset','bin2cms');
    warning on
end
if ~exist('angOffset','var')
    warning('Angle not calibrated!')
    angOffset = 0;
end
if ~exist('bin2cms','var')
    warning('Position not calibrated!')
    bin2cms = 1;
end

whl = bin2cms*whl;

dx = whl(:,1)-whl(:,3);
dy = whl(:,2)-whl(:,4);

ang = mod(atan2(dy,dx),2*pi);
da = diff(ang);
da(da<pi) = 2*pi+da(da<pi);
da(da>pi) = 2*pi - da(da>pi);
da = abs(gaussFilt(da,3,0,1));
da = gaussFilt(da,3,0,1);

dt = median(diff(t));

%dt = dt/10000;
angSpd = tsd(t,[da;0]/dt);

headVec = [dx dy];
angOffset = -angOffset;
rotMat = [[cos(angOffset);sin(angOffset)] [-sin(angOffset);cos(angOffset)]];
headVec = rotMat*headVec';
headVec = headVec';
headVec = headVec./repmat(sqrt(sum(headVec.^2,2)),[1 2]);

dxr = [diff(whl(:,1));0];
dyr = [diff(whl(:,2));0];

projSpdR = sum([dxr dyr].*headVec,2);

dxb = [diff(whl(:,3));0];
dyb = [diff(whl(:,4));0];
projSpdB = sum([dxb dyb].*headVec,2);

ix = isnan(projSpdR) & ~isnan(projSpdB);
projSpdR(ix) = projSpdB(ix);
projSpdR = abs(gaussFilt(projSpdR,3));
projSpdR = gaussFilt(projSpdR,3);

ep = ~isnan(projSpdR);
tg = t;
tg(~ep) = -1;
tg(end) = -1;
GoodRanges = thresholdIntervals(tsd(t,tg),0,'Direction','Above');

linSpd = tsd(t,projSpdR/dt);

if 0
    
   load Analysis/BehavEpochs.mat wakeEp
   linSpd = Restrict(linSpd,wakeEp);
   Xr = tsd(t,whl(:,1));
   Yr = tsd(t,whl(:,2));
   Xb = tsd(t,whl(:,3));
   Yb = tsd(t,whl(:,4));
   Xr = Data(Restrict(Xr,wakeEp));
   Yr = Data(Restrict(Yr,wakeEp));
   Xb = Data(Restrict(Xb,wakeEp));
   Yb = Data(Restrict(Yb,wakeEp));
   dx = tsd(t,headVec(:,1));
   dy = tsd(t,headVec(:,2));
   dx = Data(Restrict(dx,wakeEp));
   dy = Data(Restrict(dy,wakeEp));
   
   x = (Xr+Xb)/2;
   y = (Yr+Yb)/2;
   
   keyboard
   
   for ii=60:100:length(t)-30
       
       figure(1),clf
       plot(Xr(ii-30:ii+30),Yr(ii-30:ii+30),'ro')
       hold on
       plot(Xb(ii-30:ii+30),Yb(ii-30:ii+30),'bo')
       plot(Xr(ii-30:ii-25),Yr(ii-30:ii-25),'ko')
       plot(Xb(ii-30:ii-25),Yb(ii-30:ii-25),'ko')
       
       plot([0 dx(ii)]+x(ii),[0 dy(ii)]+y(ii),'g')
       plot(dx(ii)+x(ii),dy(ii)+y(ii),'ko')
       
       plot([0 dx(ii-30)]+x(ii-30),[0 dy(ii-30)]+y(ii-30),'r')
       plot(dx(ii-30)+x(ii-30),dy(ii-30)+y(ii-30),'ko')

       plot([0 dx(ii+30)]+x(ii+30),[0 dy(ii+30)]+y(ii+30),'b')
       plot(dx(ii+30)+x(ii+30),dy(ii+30)+y(ii+30),'ko')
       
       pause
       
       
   end
end

