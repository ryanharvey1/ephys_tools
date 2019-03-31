
% MakeWhlFile_3spots(fileBase,Session,invertSyncBool,videoRes,nSamplesPerScreen)
%
% Synchronize electrophysiological data, video streams and behavioral events
%
% Inputs
%  A .led file, a .spots file, and a .maz file
%
% Outputs
%  A .whl file and a .evt file
%
% Required parameters
%  fileBase: the (common) base name for the three input files
%
% Optional parameters
%  nSamplesPerScreen: the number of samples to display at once (default = 100)
%
% videoRes(x y): number of pixels in video input (default = [368, 240])
%
% invertSyncBool may be used to correct the polarity of the syncronization pulse
%
% Note: this is an interactive function; follow the instructions on the main window

function MakeWhlFile_3spots_GB(fileBase,invertSyncBool,videoRes,nSamplesPerScreen)
if nargin<1
    fprintf('Usage: MakeWhlFile_3spots(fileBase,Session[,invertSyncBool,videoRes,nSamplesPerScreen])\n');
    return;
end
if ~exist('nSamplesPerScreen','var') | isempty(nSamplesPerScreen),
  nSamplesPerScreen = 100;
end
if ~exist('videoRes','var') | isempty(videoRes),
  videoRes = [800,448];
end
if  ~exist('invertSyncBool','var') | isempty(invertSyncBool),
	    invertSyncBool = 0;
end

Thresh = 3; % to avoid triggers caused by 1 or 2 spurious pixels

% Load .spots, .led and .maz files

fprintf('Loading data ... ');
spots = [];
if exist([fileBase '.spots']),
  spots = load([fileBase '.spots']);
  spots0=spots; %To be used at the end
  nFrames = max(spots(:,1))+1;
end
eegAll = bload([fileBase '.led' ], [1 inf]);

if invertSyncBool %invert Sync pulses
  eeg = -eegAll(1,:);
else
  eeg = eegAll(1,:);
end
% if exist([fileBase '.maz']),
%   events = load([fileBase '.maz']);
% end
fprintf('Done\n');

size(eeg)
nFrames

% find upstrokes of eeg square wave - use mean+-1/4s.d. as trigger points
baseEeg = mean(eeg(1:2500));
gw = gausswin(1250/50);
gw = gw/sum(gw(:));
eeg = eeg(:);
eegF = convn(eeg-baseEeg,gw,'same');

syncEEG = SchmittTrigger(eegF,1000,100);
fprintf('The first EEG pulse was detected at %i sec.\n', syncEEG(1)/1250);
fprintf('The second EEG pulse was detected at %i sec.\n', syncEEG(2)/1250);


if ~isempty(spots),
  figure(1);
  clf;
  plot(spots(:,3),spots(:,4),'.','markersize',1);
  zoom on

  remove_more=1;
  while(remove_more==1)
     keyin = input('\n\n Do you want to remove some spots from analysis (yes/no)? ', 's');
     if strcmp(keyin,'yes'),
        input('In the figure window, select the area you want to delete.\n   Then click back in this window and hit ENTER...','s');
        xr = xlim;
        yr = ylim;
        [m, n] = size(spots);
 	SpotsToKeep = find(~(spots(:,3)>=xr(1) & spots(:,3)<=xr(2) & spots(:,4)>=yr(1) & spots(:,4)<=yr(2)));
 	spots = spots(SpotsToKeep,:);
	plot(spots(:,3),spots(:,4),'.','markersize',1);
      else if strcmp(keyin,'no'),
	remove_more=0;
        end
      end
  end

  fprintf('\nDetection of synchronizing LED\n------------------------------\n');
  input('   In the figure window, select the area of the sync spot so it fills the axes.\n   Then click back in this window and hit ENTER...','s');
  xr = xlim;
  yr = ylim;

  IsSyncSpot = spots(:,3)>=xr(1) & spots(:,3)<=xr(2) & spots(:,4)>=yr(1) & spots(:,4)<=yr(2);
  SyncSpots = find(IsSyncSpot); % entries in spots.txt for our LED
  PixCnt = Accumulate(spots(SyncSpots,1)+1,spots(SyncSpots,2),nFrames); % number of pixels in detected LED

  % now we get down to synchronization
  % - find flash points using "Schmitt trigger"
  % i.e. when it goes from 0 to above threshold
  syncVideo = SchmittTrigger(PixCnt,Thresh,0);

  % Before manual correction:
  %  syncEEG contains the timestamps of the SYNC upstrokes detected in the EEG
  %  syncVideo contains the timestamps of the LED onsets detected in the video
  % After manual correction:
  %  syncEEG and syncVideo will contain the corrected timestamps
  %  syncVideo will be truncated if it is longer than syncEEG, but not the other way around
  %  Therefore, length(syncVideo) <= length(syncEEG)
  fprintf('\nSynchronization of video and EEG \n--------------------------------\n');
  while 1,
      fprintf('   There are %d video flashes and %d square wave pulses.\n',length(syncVideo),length(syncEEG));
      if length(syncVideo)~=length(syncEEG),
        fprintf('   ***** MISMATCH! *****\n');
        i = input('   To manually correct this, hit ENTER; to drop flashes in excess and continue, type ''done''+ENTER - but you''d better be sure...','s');
      else
        i = input('   To manually edit this, hit ENTER; to continue, type ''done''+ENTER...','s');
      end
      if strcmp(i,'done'),
          if length(syncVideo) > length(syncEEG),
              syncVideo = syncVideo(1:length(syncEEG));
          end
          break;
      else
        [syncVideo,syncEEG] = DisplaySYNC(syncVideo,syncEEG,nSamplesPerScreen);
      end
  end
end

% Before manual correction:
%  syncEEG contains the timestamps of the SYNC upstrokes detected in the EEG
%  syncEvents contains the timestamps of the SYNC events in the .maz file
% After manual correction:
%  syncEEG and syncEvents will contain the corrected timestamps
%  syncEvents will be truncated if it is longer than syncEEG, but not the other way around
%  Therefore, length(syncEvents) <= length(syncEEG)
if ~exist('events', 'var')
    events = [];
end
if ~isempty(events),
  fprintf('\nSynchronization of events and EEG\n---------------------------------\n');
  syncEvents = events(find(events(:,2) == 83 & events(:,3) == 89));
  while 1,
      fprintf('   There are %d SYNC events and %d square wave pulses.\n',length(syncEvents),length(syncEEG));
      if length(syncEvents)~=length(syncEEG),
        fprintf('   ***** MISMATCH! *****\n');
        i = input('   To manually correct this, hit ENTER; to drop flashes in excess and continue, type ''done''+ENTER - but you''d better be sure...','s');
      else
        i = input('   To manually edit this, hit ENTER; to continue, type ''done''+ENTER...','s');
      end
      if strcmp(i,'done'),
          if length(syncEvents) > length(syncEEG),
              syncEvents = syncEvents(1:length(syncEEG));
          end
          break;
      else
        [syncEvents,syncEEG] = DisplaySYNC(syncEvents,syncEEG,nSamplesPerScreen);
      end
  end
end

if ~isempty(spots),
  figure(1);
  subplot(2,2,1);
  plot(diff(syncVideo),'.-')
  ylabel('# video frames between flashes');
  xlabel('flash #');

  subplot(2,2,2)
  plot(diff(syncEEG), '.-')
  ylabel('# EEG samples between flashes');
  xlabel('flash #');

  subplot(2,2,3);
  [b bint r] = regress(syncEEG(1:length(syncVideo)),[syncVideo,ones(size(syncVideo))]);
  plot(r/1.25,'.')
  xlabel('Flash #');
  ylabel('deviation from linear fit (ms)');

  subplot(2,2,4);
  hold off;
  plot(diff(syncVideo)./diff(syncEEG(1:length(syncVideo)))*1250,'.');
  FilterLen = 10;
  f = filter(ones(FilterLen,1)/FilterLen, 1,diff(syncVideo)./diff(syncEEG(1:length(syncVideo)))*1250);
  hold on;
  plot(FilterLen:length(f),f(FilterLen:end),'r');
  ylabel('Frame rate (red is smoothed)');
  xlabel('Flash #');

  % now align them - any outside sync range become NaN
  FrameSamples = interp1(syncVideo,syncEEG(1:length(syncVideo)),1:nFrames,'linear',NaN);

  % THE FOLLOWING ONLY WORKS IF YOU HAVE 2 LEDS

  % find non-sync spots
  NonSyncSpots = find(~IsSyncSpot);
  % find those with 2 non-sync spots per frame
  NonSyncCnt = Accumulate(1+spots(NonSyncSpots,1),1);
  GoodSpots = NonSyncSpots(find(NonSyncCnt(1+spots(NonSyncSpots,1))==2));
  % select front versus rear LED in color space
  figure(2);
  clf;
  rgb = ycbcr2rgb(spots(GoodSpots,7:9)/256);
  plot(rgb(:,1),rgb(:,3),'.','markersize',1);
  xlabel('Red');
  ylabel('Blue');
  xlim([-0.1 1.1])
  ylim([-0.1 1.1])
  %zoom on;
  
  fprintf('\nDiscrimination of front and rear lights\n---------------------------------------\n');
  fprintf('\n Select the front spot');
  %input('   In the figure window, select the front spot so it fills the axes.\n   Then click back in this window and hit ENTER...','s');
  %xr = xlim;
  %yr = ylim;
  %IsFrontSpot = rgb(:,2)>=xr(1) & rgb(:,2)<=xr(2) & rgb(:,3)>=yr(1) & rgb(:,3)<=yr(2);
  IsFrontSpot = ClusterPoints(rgb(:,[1 3]),0);
  FrontSpots = GoodSpots(find(IsFrontSpot));
  RearSpots = GoodSpots(find(~IsFrontSpot));
  
  figure(2);
  subplot(2,1,1);hold on;
  plot(spots(FrontSpots,3),spots(FrontSpots,4),'.','color',[1 0 0],'markersize',10,'linestyle','none');
  plot(spots(RearSpots,3),spots(RearSpots,4),'.','color',[0 1 0],'markersize',10,'linestyle','none');
  subplot(2,1,2);hold on;
  k=0;

	% now make trajectory
  HeadPos = -1*ones(nFrames,4);
  HeadPos(1+spots(FrontSpots,1),1:2) = spots(FrontSpots,3:4);
  HeadPos(1+spots(RearSpots,1),3:4) = spots(RearSpots,3:4);
  % interpolate missing stretches up to 10 frames long
  cHeadPos = CleanWhl(HeadPos, 10, inf);
  cHeadPos(find(cHeadPos==-1))=NaN; % so it doesn't interpolate between 100 and -1 and get 50.

	% now make wheel file by interpolating
  TargetSamples = 0:32:length(eeg);
  GoodFrames = find(isfinite(FrameSamples));
  Whl(:,1:4) = interp1(FrameSamples(GoodFrames),cHeadPos(GoodFrames,:),TargetSamples,'linear',-1);
  Whl(find(~isfinite(Whl)))=-1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   WhlRL = [Whl,zeros(size(Whl,1),1)];
%   WhlSessioNum = zeros(size(Whl,1),1);
% 
%   N_session = size(Session,1);
% 
%   for ii=1:1:N_session
%      SessionStart = Session(ii,1);
%      SessionEnd   = Session(ii,2);
%      SessionRL    = Session(ii,3);
% 
%      SessionStart = floor(SessionStart * 39.0625);
%      SessionEnd = floor(SessionEnd * 39.0625);
% 
%      WhlRL(SessionStart:SessionEnd,5)=SessionRL;
%      WhlSessionNum(SessionStart:SessionEnd,1)=ii;
%   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
size(Whl)
k=0;
plotWhl = Whl(find(Whl(:,1)~=-1),:);
fprintf('\nThe Whl file data\n----------------------\n');
%while ~strcmp(input('   Hit ENTER to show the next 100 samples, or type ''done''+ENTER to proceed...','s'),'done'),
%    k = k+1;
%    if k*100 > length(plotWhl), break; end
%    cla;
%    set(gca,'xlim',[0 videoRes(1)],'ylim',[0 videoRes(2)]);
%    plot(plotWhl((k-1)*100+1:k*100,1),plotWhl((k-1)*100+1:k*100,2),'.','color',[1 0 0],'markersize',10,'linestyle','none');
%    plot(plotWhl((k-1)*100+1:k*100,3),plotWhl((k-1)*100+1:k*100,4),'.','color',[0 1 0],'markersize',10,'linestyle','none');
%    for j=(k-1)*100+1:k*100, line([plotWhl(j,1) plotWhl(j,3)],[plotWhl(j,2) plotWhl(j,4)],'color',[0 0 0]); end
%    set(gca,'xlim',[0 videoRes(1)],'ylim',[0 videoRes(2)]);
%end

% while ~strcmp(input('   Hit ENTER to show the next session, or type ''done''+ENTER to proceed...','s'),'done'),
%     k = k+1;
%     if k > N_session, break; end
%     cla;
%     set(gca,'xlim',[0 videoRes(1)],'ylim',[0 videoRes(2)]);
%     findPlotSession=find(WhlSessionNum(:,1)==k);
%     %findPlotSessionGood=find(Whl(findPlotSession,1)~=-1);
%     plot(Whl(findPlotSession,1),Whl(findPlotSession,2),'.','color',[1 0 0],'markersize',10,'linestyle','none');
%     plot(Whl(findPlotSession,3),Whl(findPlotSession,4),'.','color',[0 1 0],'markersize',10,'linestyle','none');
%     %for j=1:length(findPlotSessionGood)
%     %    line([Whl(findPlotSessionGood(j),1) Whl(findPlotSessionGood(j),3)],[Whl(findPlotSessionGood(j),2) Whl(findPlotSessionGood(j),4)],'color',[0 0 0]); 
%     %end
%     set(gca,'xlim',[0 videoRes(1)],'ylim',[0 videoRes(2)]);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~isempty(events),
  % Update events: remove all initial SYNC events and replace them with the corrected ones (according to user input)
  nonSYNC = events(events(:,2) ~= 83 | events(:,3) ~= 89,:);
  events = [nonSYNC;[syncEvents 83*ones(size(syncEvents)) 89*ones(size(syncEvents))]];
  events = sortrows(events,1);
  % Throw any events occurring after the last SYNC event
  lastSYNC = find(events(:,2) == 83 & events(:,3) == 89);lastSYNC = lastSYNC(end);
  events = events(1:lastSYNC,:);
  % Synchronize events on electrophysiological data
  timestamps = interp1(syncEvents,syncEEG(1:length(syncEvents))/1250*1000,events(:,1),'linear',-1);
  events(:,1) = timestamps;
  events = double(uint32(events));
end

if ~isempty(spots),
  figure(1);
end

figure(2)
  clf;
  plot(spots(:,3),spots(:,4),'.','markersize',1);

done=0;
while ~done
    fprintf('Draw a line on the figure\n')
    gg = ginput(2);
    lpixel = norm(gg(1,:)-gg(2,:));
    lcms = input('What is the actual size of this line (in cms)','s');
    ratio = str2num(lcms)/lpixel;
    fprintf('Ratio cms/pixel is: %f\n',ratio)
    answer = input('Satisfied? [Y/N]','s');
    keyboard
    if strcmp(answer,'Y')
        done=1;
    elseif strcmp(answer,'N')
        fprintf('OK, then let''s restart!\n')
    else
        fprintf('Don''t understand the answer, so we restart\n')
    end
end

Whl = ratio*Whl;

while 1,
  i = input('\n\nSave to disk (yes/no)? ', 's');
  if strcmp(i,'yes') | strcmp(i,'no'), break; end
end
if i(1) == 'y'
  if ~isempty(spots),
    fprintf('Saving %s\n', [fileBase '.whl']);
    msave([fileBase '.whl'], Whl);
%     fprintf('Saving %s\n', [fileBase '.whlrl']);
%     msave([fileBase '.whlrl'], WhlRL);
  end
  if ~isempty(events),
    fprintf('Saving %s\n', [fileBase '.evt']);
    msave([fileBase '.evt'],events);
  end
end

if 0 % no longer needed because neuroscope is cool
  while 1,
  i = input('\n\nSave all position information as jpeg (yes/no)? ', 's');
  if strcmp(i,'yes') | strcmp(i,'no'), break; end
  end
  if i(1)=='y'
    figure;
    plot(plotWhl(:,1),plotWhl(:,2), '.','color',[0.5 0.5 0.5],'markersize',10,'linestyle','none' );
    set(gca, 'xlim', [0 videoRes(1)], 'ylim', [0 videoRes(2)], 'Position', [0 0 1 1]);
    set(gca, 'color', [0 0 0]);
    set(gcf, 'inverthardcopy', 'off')
    print(gcf, '-djpeg100', [fileBase '.jpg']);
    eval(['!convert -geometry ' num2str(videoRes(1)) 'x' num2str(videoRes(2)) ' ' fileBase '.jpg' ' ' fileBase '.jpg']);
  end
  end
