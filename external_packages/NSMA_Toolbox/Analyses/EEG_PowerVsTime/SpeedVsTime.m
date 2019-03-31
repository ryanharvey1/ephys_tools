function SpeedVsTime(X,Y,epochs,varargin)

% SpeedVsTime  Plots animal's (smoothed) speed over time, based on tracker data
%
% function SpeedVsTime(X,Y,epochs,varargin)
%
%  X,Y = position tsd's (pixels)
%  epochs = struct with names and time intervals of epochs
%  varargin (optional):  xlims = x-axis limits -- 1x2 array of start & end timestamps
%
%  Produces plot of speed (pixels/sec) vs. time, with epoch boundaries shown
%
%  PL '02, last modified '04 by MN


if ~isempty(varargin)
    xlims = varargin{1};
else
    xlims = [epochs.intervals{1}(1) epochs.intervals{end}(2)];
end %if

xsec = range(X,'sec');
xdata = data(X);
ydata = data(Y);
dt = diff(xsec);

% take out dt zeros
xsec(find(dt == 0)+1,:)=[]; % remove 2nd frame of diff pair
xdata(find(dt == 0)+1,:)=[];
ydata(find(dt == 0)+1,:)=[];
X = tsd(xsec*10000,xdata);
Y = tsd(xsec*10000,ydata);


% get cumulative distance travelled
dx = diff(xdata);
dy = diff(ydata);
dt = diff(xsec);
distances = sqrt(dx.^2 + dy.^2);
cumdist = [0; cumsum(distances)];
cumdist_tsd = tsd([0; dt*10000], cumdist);


% "pad" (add 4 sec to) beginning and end of the data so that the smoother function below doesn't
% create strange results at beginning or end of the data
ExtraTimesBegin = linspace(xsec(1)-4,xsec(1),241)';
ExtraTimesBegin = ExtraTimesBegin(1:end-1);
ExtraDataBegin = ones(240,1)*cumdist(1);

ExtraTimesEnd = linspace(xsec(end),xsec(end)+4,241)';
ExtraTimesEnd = ExtraTimesEnd(2:end);
ExtraDataEnd = ones(240,1)*cumdist(end);

xsec = [ExtraTimesBegin; xsec; ExtraTimesEnd];
cumdist = [ExtraDataBegin; cumdist; ExtraDataEnd];
dt=diff(xsec);


% pad timegaps to eliminate strange effects in subsequent smoothing
% flag gaps in the time sequence longer than 2 sec
igaps=find(dt > 2); 
    
if ~isempty(igaps)

    xsec_before_gaps = zeros(length(xsec),1);  
    xsec_before_gaps(igaps) = 1;  
    xsec_after_gaps = zeros(length(xsec),1);  
    xsec_after_gaps(igaps+1) = 1;  

    % add 60 frames after frame before gap, and 60 frames before frame after gap
    for i = 1:60
        ExtraTimesAfterStops(:,i) = xsec(find(xsec_before_gaps == 1))+(i/60);
        ExtraTimesBeforeStarts(:,61-i) = xsec(find(xsec_after_gaps == 1))-(i/60); 
    end
    ExtraTimesAfterStops = reshape(ExtraTimesAfterStops',length(igaps)*60,1);
    ExtraTimesBeforeStarts = reshape(ExtraTimesBeforeStarts',length(igaps)*60,1);

    ExtraDataAfterStops = repmat(cumdist(find(xsec_before_gaps == 1)),1,60);
    ExtraDataBeforeStarts = repmat(cumdist(find(xsec_after_gaps == 1)),1,60);

    ExtraDataAfterStops = reshape(ExtraDataAfterStops',length(igaps)*60,1);
    ExtraDataBeforeStarts = reshape(ExtraDataBeforeStarts',length(igaps)*60,1);

    % may bomb on these big matrices
    mat = [[xsec; ExtraTimesAfterStops; ExtraTimesBeforeStarts] [cumdist; ExtraDataAfterStops; ExtraDataBeforeStarts] ...
            [zeros(length(xsec),1); ones(length([ExtraTimesAfterStops; ExtraTimesBeforeStarts]),1)]];
    mat = sortrows(mat);
    
else  
    mat = [xsec cumdist];
    
end % if ~isempty...


% Smooth cumulative position data (2 sec hamming window) in deg/frame (60 frames/sec)
% (only roughly 2 sec - more accurate just to say "120-frame" hamming win. - esp. when
% there are timegaps in the data)
cumdist_tsd=tsd(mat(:,1)*10000,mat(:,2));
Xs = SmoothTSD(cumdist_tsd,120); clear cumdist_tsd;


% get rid of padded data
if ~isempty(igaps)
    Xs = [range(Xs,'ts') data(Xs) mat(:,3)];
else
    Xs = [range(Xs,'ts') data(Xs)];
end
% get rid of data added to beginning and end of session
Xs = Xs(241:end-240,:);
% get rid of data added in gaps
if ~isempty(igaps)
    Xs(find(Xs(:,3) == 1),:) = [];
end


% calculate smoothed speed (pixels per second)
diff_sec=diff(Xs(:,1)/10000);
diff_data=diff(Xs(:,2));
speed=[0; diff_data./diff_sec];
sec = Xs(:,1)/10000;


% flag gaps in the time sequence longer than 2 sec
igaps=find(diff_sec > 2);
unknown_data=ones(length(speed),1)*NaN;
unknown_data(igaps)=0;
unknown_data(igaps+1)=0;
speed(igaps+1)=NaN;


%PLOT

%plot speed vs. time
H = plot(sec/60, speed, '-');

%mark gaps in data with cyan blue
hold on; plot(sec/60, unknown_data, 'c');
plot([xlims(1)/600000; sec(1)/60],[0;0],'c');
plot([sec(end)/60; xlims(end)/600000],[0;0],'c');

%plot epoch lines
hold on;
t0 = xlims(1)/600000;
t1 = xlims(2)/600000;
%y1=max(speed);
y1=prctile(speed,99.9);
%y0 = 0;
y0 = -y1/10;
DrawEpochLines(y1,y0,epochs);
    
%axes
axis([t0 t1 y0 y1]);
set(H(1), 'Color', [.75 .75 .75]);
set(H(1), 'LineWidth', 0.1);
str(1) = {'speed'};
str(2) = {'(pixels/sec)'};
h = ylabel(str,'Units','normalized');
g=get(h,'Position');
g(1) = g(1)-.025;
ylabel(str,'Position',g,'Units','normalized');


%--------------------------------------------
function DrawEpochLines(y0,y1,epochs)

    for iep = 1:length(epochs.names)
        t0 = epochs.intervals{iep}(1)/10000/60;
        t1 = epochs.intervals{iep}(2)/10000/60;
        tm = (t0+t1)/2;
        hhl = line([t0 t0], [y0 y1]);
        set(hhl,'Color','r');
        hhl = line([t1 t1], [y0 y1]);
        set(hhl,'Color','b');
        text(tm,y0-(y1-y0)/15,epochs.names{iep});
    end