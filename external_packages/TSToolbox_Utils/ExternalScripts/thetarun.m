%function thrun = thetarun(FileBase,Periods,frin,frout,ch2use)
% runs automatic Gaussian HMM segmentation of the period of recording
% given by Periods in seconds
function thrun = thetarun(FileBase,varargin)
load([FileBase '.eegseg.mat']);
Par = LoadXml([FileBase '.xml']);
[Periods,frin,frout,ch2use] = DefaultArgs(varargin,{ [t(1) t(end)]*Par.lfpSampleRate, [5 11], [1 5; 11 15],1});
%%%%%%%%%%%%%%%%%%%%
%PARAMS - temporary
%per2use =1;
%ch2use = 1;
%%%%%%%%%%%%%%%%%%

%frout = frout(per2use,:);

%find freq. indexes within and outsize the band of interest
thfin = findind(f,frin);
thfout = findind(f,frout);

timein = findind(t,Periods);
timeout =setdiff([1:length(t)],timein);
thratio = log(mean(squeeze(y(timein,thfin,ch2use)),2))-log(mean(squeeze(y(timein,thfout,ch2use)),2));
nStates =2;
% fit gaussian mixture and to HMM - experimental version .. uses only
% thetaratio 
[theState thhmm thdec] = gausshmm(thratio,nStates,1,0);

%find which state is theta  and make step vector from 0 to 1 (outside to
%inside)
for ii=1:nStates 
	thratio_st(ii) = mean(thratio(find(theState==ii)));
end
[dummy TheInd] = max(thratio_st);
InTh = zeros(length(t),1);
InTh(timein) = (theState==TheInd);
InTh(timeout) = 0; % that puts all times outside the periods as not theta

%detect beg and end of the theta runs
thebeg= SchmittTrigger(InTh,0.9, 0.9);
theend= SchmittTrigger(-InTh,-0.9, -0.9);

if thebeg(1)>theend(1)
    theend =theend(2:end);
end

if thebeg(end)>theend(end)
    thebeg =thebeg(1:end-1);
end
%theend = theend-1; %shif tend back as the it is detected on donw/up

thrun =round([t(thebeg) t(theend)]*Par.lfpSampleRate);

if nargout<1
    figure
    myf= [min([frin(:);frout(:)]) max([frin(:);frout(:)])];
    myfi = find(f>myf(1) & f<myf(2));
    imagesc(t,f,log(sq(y(:,myfi,1)))');
    axis xy
    ca = caxis; caxis([ca(1) 0.8*ca(2)]);
    hold on
    h1 = Lines(thrun(:,1)/1250,[],'b');
    h2 = Lines(thrun(:,2)/1250,[],'r');
    h3 = Lines(Periods(:),[],'k');
	set(h1,'LineWidth',2);
	set(h2,'LineWidth',2);
    set(h3,'LineWidth',2);
    TimeBrowse(100,100);
 %   msave([FileBase '.thrun'],'thrun');
end

% nSeg = size(thebeg,1);
% for ii=1:nSeg
%     thsegm(ii) = mean(thratio([thebeg(ii):theend(ii)]));
%     thsegv(ii) = std(thratio([thebeg(ii):theend(ii)]));
% end
%     
% pow = log(sq(y(:,myfi,1)));
% mpow = mean(pow(find(InTh),:));
% for ii=1:nSeg
%     md = mahal(mpow,pow([thebeg(ii):theend(ii)],:));
%     
%     
%     
    

function ind = findind(f,per)
ind =[];f=f(:);
for ii=1:size(per,1)
   ind = [ind; find(f>per(ii,1) & f<per(ii,2))];
end
return