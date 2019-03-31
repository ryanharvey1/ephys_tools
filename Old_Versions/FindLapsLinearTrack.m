function [laps,right,left]=FindLapsLinearTrack(data_video,data_video_spk)
% FindLapsLinearTrack: This function locates and splits up laps on the linear track
% 
% Input: 
%           data_video: tracking data matrix containing the following structure*
% 
%           data_video = [timestamp,x,y,...;
%                         timestamp,x,y,...;
%                         timestamp,x,y,...;
%                         .................]
% 
%           * note that if data_video is not sampled at 30hz, change the
%           smoothing parameters
% 
% Output:
%           laps: cell array containing each individual lap seperated from data_video.
%                 Each row is a lap & each column is a direction (column 1 right, column 2 left)
%           right: right laps
%           left: left laps
% 
% Ryan Harvey 2018
% 
% 
% SPLIT UP TRACK BY LAP
% smooth x coordinate 
% x=smooth(data_video(:,2),450);

x=smoothdata(data_video(:,2));



% locate angle of lap to pull apart right and left laps
A=XYangle(rescale(data_video(:,1),0,1),rescale(x,0,1));

% locate initial starts and end points by seperating path by above angle
rightidx=A>200;
[startr,endsr,ngroupsr]=findgroups(rightidx);
leftidx=A<200;
[startl,endsl,ngroupsl]=findgroups(leftidx);

% a full lap must traverse the entire track
% locate lower and higher cut offs
xrange=range([min(data_video(:,2)),max(data_video(:,2))])*.35;
lowercut=min(data_video(:,2))+xrange;
highercut=max(data_video(:,2))-xrange;

% find laps that traverse the entire track by the above criteria
while true % repeat until lap number stops decreasing
    % if lap does not cross middle, assign that lap to the opposite direction
    % right
    for i=1:length(startr)
        if (any(data_video(startr(i):endsr(i),2)>highercut) && any(data_video(startr(i):endsr(i),2)<lowercut))==0
            leftidx(startr(i):endsr(i))=ones(1,length(startr(i):endsr(i)));
            rightidx(startr(i):endsr(i))=zeros(1,length(startr(i):endsr(i)));
        end
    end
    % left
    for i=1:length(startl)
        if (any(data_video(startl(i):endsl(i),2)>highercut) && any(data_video(startl(i):endsl(i),2)<lowercut))==0
            rightidx(startl(i):endsl(i))=ones(1,length(startl(i):endsl(i)));
            leftidx(startl(i):endsl(i))=zeros(1,length(startl(i):endsl(i)));
        end
    end
    
    ngroupsr_old=ngroupsr;
    ngroupsl_old=ngroupsl;
    [startr,endsr,ngroupsr]=findgroups(rightidx);
    [startl,endsl,ngroupsl]=findgroups(leftidx);
    
    % when the number of laps stops decreasing, break out of while loop
    if ngroupsr==ngroupsr_old && ngroupsl==ngroupsl_old
        break
    end
end

% Fix the start and stop points of each lap. This is done because
% calculating the angle on smoothed data will not accurately locate the
% exact peaks and valley index

% fix points after peaks
% toend=endsr(end)>endsl(end);
for i=1:length(startl)-1
    
    lpI=length(data_video(startl(i):endsl(i),2))+(startl(i)-1);
    [~,mI]=max(data_video(startl(i):endsl(i),2));
    mI=mI+(startl(i)-1);
    
    if mI<lpI
        leftidx(mI+1:lpI)=zeros(1,length(mI+1:lpI));
        rightidx(mI+1:lpI)=ones(1,length(mI+1:lpI));
    end
end
[startr,endsr,~]=findgroups(rightidx);
[startl,endsl,~]=findgroups(leftidx);

% fix points before valley
for i=2:length(startl)
    
    lpI=startl(i);
    [~,mI]=min(data_video(startl(i):endsl(i),2));
    mI=mI+(startl(i)-1);
    
    if lpI<mI
        leftidx(lpI:mI+1)=zeros(1,length(lpI:mI+1));
        rightidx(lpI:mI+1)=ones(1,length(lpI:mI+1));
    end
end
[startr,endsr,~]=findgroups(rightidx);
[startl,endsl,~]=findgroups(leftidx);

% fix points after valley
toend=endsr(end)>endsl(end);
for i=1:length(startr)-toend

    lpI=endsr(i);
    [~,mI]=min(data_video(startr(i):endsr(i),2));
    mI=mI+(startr(i)-1);
    
    if mI<lpI
        leftidx(mI+1:lpI)=ones(1,length(mI+1:lpI));
        rightidx(mI+1:lpI)=zeros(1,length(mI+1:lpI));
    end
end
[startr,endsr,~]=findgroups(rightidx);
[startl,endsl,~]=findgroups(leftidx);

% fix points before peak
for i=2:length(startr)
    
    lpI=startr(i);
    [~,mI]=max(data_video(startr(i):endsr(i),2));
    mI=mI+(startr(i)-1);
    
     if lpI<mI
        leftidx(lpI:mI+1)=ones(1,length(lpI:mI+1));
        rightidx(lpI:mI+1)=zeros(1,length(lpI:mI+1));
    end
end
[startr,endsr,~]=findgroups(rightidx);
[startl,endsl,groupsl]=findgroups(leftidx);

figure
plot(data_video(rightidx,1),data_video(rightidx,2),'.g');hold on
plot(data_video(leftidx,1),data_video(leftidx,2),'.r')
title('Cleaned lap split')


for i=1:length(startl)
    c=corr(data_video(startl(i):endsl(i),1),data_video(startl(i):endsl(i),2))
end

figure
for i=1:length(startr)
    c=corr(data_video(startr(i):endsr(i),1),data_video(startr(i):endsr(i),2))
     plot(data_video(startr(i):endsr(i),1),data_video(startr(i):endsr(i),2));hold on
    if c>0
        test=1
    end
end

% store data in a cell array and in lap vectors: right and left
right=[];
for i=1:length(startr)
    laps{i,1}=data_video_spk(find(data_video_spk(:,1)==data_video(startr(i),1)):...
        find(data_video_spk(:,1)==data_video(endsr(i),1)),:);
    right=[right;data_video_spk(find(data_video_spk(:,1)==data_video(startr(i),1)):...
        find(data_video_spk(:,1)==data_video(endsr(i),1)),:)];
end

left=[];
for i=1:length(startl)
    laps{i,2}=data_video_spk(find(data_video_spk(:,1)==data_video(startl(i),1)):...
        find(data_video_spk(:,1)==data_video(endsl(i),1)),:);
    left=[left;data_video_spk(find(data_video_spk(:,1)==data_video(startl(i),1)):...
        find(data_video_spk(:,1)==data_video(endsl(i),1)),:)];
end
end

% /////////// Old code //////////////

% % plot feature for debugging
% plots=1;
% 
% % remove any nans
% data_video_spk(isnan(data_video_spk(:,1)),:)=[];
% 
% % identify x
% x=data_video_spk(:,2);
% 
% % find peak above 85% of max peak with min distance between peaks at least 11 seconds
% [pktp,lctp]=findpeaks(x,'MinPeakHeight',max(x)*.85,'MinPeakDistance',330);
% 
% % find valley below 15% of max valley with min distance between valleys at least 11 seconds
% [pktv,lctv]=findpeaks(-x,'MinPeakHeight',-(max(x)*.15),'MinPeakDistance',330);
% 
% if plots==1;figure;plot(x,'.k');hold on;scatter(lctp,pktp,'r');scatter(lctv,-pktv,'g');hold off;end
% 
% % if first peak is before the first valley, flip the path so that valleys are first
% if lctp(1)<lctv(1)
% %     x=rescale(-x,0,1);
% x=-x;
%     [pktp,lctp]=findpeaks(x,'MinPeakHeight',max(x)*.85,'MinPeakDistance',330);
%     [pktv,lctv]=findpeaks(x,'MinPeakHeight',(min(x)*.85),'MinPeakDistance',330);
%     
%     if plots==1;figure;plot(x,'.k');hold on;scatter(lctp,pktp,'r');scatter(lctv,-pktv,'g');hold off;end
% end
% 
% % check to find multiple valleys exist before first peak
% if sum(lctv<lctp(1))>1
%     [~,I]=min(-pktv(lctv<lctp(1)));
%     lctv=lctv(I:end);
%     pktv=pktv(I:end);
% end
% 
% [m,I]=max([lctv(end),lctp(end)]);
% if I==1
%     lctv(end)=length(data_video_spk)
% elseif I==2
%     lctp(end)=length(data_video_spk)
% end
% 
% 
% % first left lap
% lap{1,1}=data_video_spk(1:lctp(1),:);
% 
% if plots==1;figure;plot(data_video_spk(1:lctp(1),1),data_video_spk(1:lctp(1),2),'.r');hold on;end
% 
% % Left
% i=1;
% while true
%     if i+1>length(lctp)
%         break
%     end
%     lap{i+1,1}=data_video_spk(lctv(i+1):lctp(i+1),:);        
%     if plots==1;plot(data_video_spk(lctv(i+1):lctp(i+1),1),data_video_spk(lctv(i+1):lctp(i+1),2),'.r');hold on;end
%     i=i+1;
% end
% 
% % Right
% i=1;
% while true
%     if i>length(lctp)
%         break
%     end
%     lap{i,2}=data_video_spk(lctp(i):lctv(i+1),:);
%     if plots==1;plot(data_video_spk(lctp(i):lctv(i+1),1),data_video_spk(lctp(i):lctv(i+1),2),'.g');hold on;end
%     i=i+1;
% end
% 
% Left=[];
% for i=1:size(lap,1)
%     Left=[Left;lap{i,1}];
% end
% 
% 
% Right=[];
% for i=1:size(lap,1)
%     Right=[Right;lap{i,2}];
% end
% 
% close all




% 
% clear lapL lapR
% if lctp(1)<lctv(1)
%     lap{1}=data_video_spk(1:lctv(1),:);
% elseif lctp(1)>lctv(1)
%     lapL{1}=data_video_spk(1:lctp(1),:);
% %     plot(data_video_spk(1:lctp(1),1),data_video_spk(1:lctp(1),2),'.g');hold on
% 
%     if sum(lctv<lctp(1))>1
%         [~,I]=min(-pktv(lctv<lctp(1)));
%         lctv=lctv(I:end);
%         pktv=pktv(I:end);
%     end
%     i=2;
%     while true
%         try
%             lapL{i}=data_video_spk(lctp(i-1):lctv(i),:);
% %             plot(data_video_spk(lctp(i-1):lctv(i),1),data_video_spk(lctp(i-1):lctv(i),2),'.r');hold on
%             
%             lapR{i-1}=data_video_spk(lctv(i):lctp(i+1),:);
% %             plot(data_video_spk(lctv(i):lctp(i+1),1),data_video_spk(lctv(i):lctp(i+1),2),'.g');hold on
%             i=i+1;
%         catch
%             break
%         end
%     end
% end
% figure
% for i=1:length(lapR)
%     plot(lapR{i}(:,1),lapR{i}(:,2),'.g');hold on
% end
% for i=1:length(lapL)
%     plot(lapL{i}(:,1),lapL{i}(:,2),'.r');hold on
% end
% 
% 
% 
% 
% %        findpeaks(x(start:ends),'MinPeakHeight',.85)
% %%
% figure;
% x=rescale(data_video_spk(:,2),0,1);
% 
% [start,ends,ngroups]=findgroups(contiguousframes(x>.85,12))
% 
% 
% % find peak
% start=[];
% ends=[];
% npeak=1;
% tspeak=[];
% for i=1:length(x)
%     if x(i)>.85 && isempty(start)
%         start=i;
%     end
%     if x(i)<.85 && ~isempty(start)
%         ends=i;
%     end
%     if ~isempty(start) && ~isempty(ends)
%         plot(data_video_spk(start:ends,1),x(start:ends));hold on
%         tspeak=[tspeak;data_video_spk(start:ends,1)];
%         peak(npeak)=max(x(start:ends));
%         npeak=npeak+1;
%         start=[];
%         ends=[];
%     end
% end
% 
% split=tspeak(find(diff(tspeak)>.5e7));
% 
% for i=1:length(split)
%     plot(x(1:find(split(i)==data_video_spk(:,1))))
%     if i==1
%         [m(i),I(i)]=max(x(1:find(split(i)==data_video_spk(:,1))));
%     else
%         [m(i),I(i)]=max(x(find(split(i-1)==data_video_spk(:,1)):find(split(i)==data_video_spk(:,1))));
%     end
% end
% 
% peakrow=find(x==m);
% intersect(x,m)
% scatter(data_video_spk(peakrow,1),x(peakrow))
% 
% 
% % find valley
% start=[];
% ends=[];
% nvalley=1;
% for i=1:length(x)
%     if x(i)>.15 && isempty(start)
%         start=i;
%     end
%     if x(i)<.15 && ~isempty(start)
%         ends=i;
%     end
%     if ~isempty(start) && ~isempty(ends)
%         plot(data_video_spk(start:ends,1),x(start:ends));hold on
%         tsvalley{nvalley}=data_video_spk(start:ends,1);
%         valley(nvalley)=min(x(start:ends));
%         nvalley=nvalley+1;
%         start=[];
%         ends=[];
%     end
% end
