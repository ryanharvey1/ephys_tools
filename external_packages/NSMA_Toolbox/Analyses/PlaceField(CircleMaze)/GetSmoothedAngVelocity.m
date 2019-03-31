function [dps_tsd] = GetSmoothedAngVelocity(V, horseshoe_flag, varargin);

% GetSmoothedAngVelocity  Creates tsd object of smoothed angular velocity data (degrees per second)
%
% [dps_tsd] = GetSmoothedAngVelocity(V, horseshoe_flag, varargin);
%
% INPUTS:
%   V - tsd of position data on circular track (degrees)
%   horseshoe_flag - 1=animal runs in "horseshoe" pattern on circular track, where it runs once
%       around circle in one direction, turns and runs once around in reverse direction, etc.
%       0=animal runs regular laps around circle
%   VARARGIN parameter:
%       if "1", V is tsd in form of CUMULATIVE degrees rat has traveled.  Data is
%       in form [0; cumulative diffs].
% OUTPUTS:
%   dps_tsd - tsd of smoothed angular velocity data (degrees/sec)
%
% MN 6/03


% take out dt zeros
Vsec = Range(V,'sec');
Vdata = data(V);
dt=diff(Vsec);
Vsec(find(dt == 0)+1,:)=[]; % remove 2nd frame of diff pair
Vdata(find(dt == 0)+1,:)=[];
%V = tsd(Vsec*10000,Vdata);

if length(varargin) == 0
% GET CUMULATIVE POSITION DATA (cumulative degrees)
    dt=diff(Vsec); 
    dV=diff(Vdata);
    ix_crossings = find(abs(dV) > 300);

    if horseshoe_flag == 1  
        % correct dV when rat reverses direction at transition between end of lap and beginning of next
        DistFromLapEnd_before = abs(Vdata(ix_crossings)-360);
        DistFromLapEnd_after = abs(Vdata(ix_crossings+1)-(-360));
        corrected_dV_vals = DistFromLapEnd_after-DistFromLapEnd_before;
        dV(ix_crossings) = corrected_dV_vals;

        % correct dV when rat reverses direction at middle of lap
        markhalflaps = ones(size(Vdata));
        markhalflaps(find(Vdata < 0)) = 0; %flag data in 1st 1/2-lap with "0" & in 2nd 1/2-lap with "1"
        markhalflapsplus1 = [markhalflaps(2:end); 0];
        ixbeforemiddle = find(markhalflaps == 0 & markhalflapsplus1 == 1);
        DistFromMiddle_before = abs(Vdata(ixbeforemiddle));
        DistFromMiddle_after = Vdata(ixbeforemiddle+1);
        corrected_dV_vals = DistFromMiddle_after-DistFromMiddle_before;
        dV(ixbeforemiddle) = corrected_dV_vals;
    else
        bad_dV_vals = dV(ix_crossings);
        corrected_dV_vals = -sign(bad_dV_vals).*(360-abs(bad_dV_vals));
        dV(ix_crossings) = corrected_dV_vals;
    end %if horseshoe

    cvs = [0; cumsum(dV)];
end %if length(varargin)

if ~isempty(varargin) 
    if varargin{1} == 1
        cvs = data(V);
        Vsec = range(V,'sec');
    end %if varargin
end %if ~isempty


% "pad" (add 4 sec to) beginning and end of the data so that the smoother function below doesn't
% create strange results at beginning or end of the data
ExtraTimesBegin = linspace(Vsec(1)-4,Vsec(1),241)';
ExtraTimesBegin = ExtraTimesBegin(1:end-1);
ExtraDataBegin = ones(240,1)*cvs(1);

ExtraTimesEnd = linspace(Vsec(end),Vsec(end)+4,241)';
ExtraTimesEnd = ExtraTimesEnd(2:end);
ExtraDataEnd = ones(240,1)*cvs(end);

Vsec = [ExtraTimesBegin; Vsec; ExtraTimesEnd];
cvs = [ExtraDataBegin; cvs; ExtraDataEnd];
dt=diff(Vsec);

% pad timegaps to eliminate strange effects in subsequent smoothing
% flag gaps in the time sequence longer than 2 sec
igaps=find(dt > 2);
    
if ~isempty(igaps)

    Vsec_before_gaps = zeros(length(Vsec),1);  
    Vsec_before_gaps(igaps) = 1;  
    Vsec_after_gaps = zeros(length(Vsec),1);  
    Vsec_after_gaps(igaps+1) = 1;  

    % add 60 frames after frame before gap, and 60 frames before frame after gap
    for i = 1:60
        ExtraTimesAfterStops(:,i) = Vsec(find(Vsec_before_gaps == 1))+(i/60);
        ExtraTimesBeforeStarts(:,61-i) = Vsec(find(Vsec_after_gaps == 1))-(i/60); 
    end
    ExtraTimesAfterStops = reshape(ExtraTimesAfterStops',length(igaps)*60,1);
    ExtraTimesBeforeStarts = reshape(ExtraTimesBeforeStarts',length(igaps)*60,1);

    ExtraDataAfterStops = repmat(cvs(find(Vsec_before_gaps == 1)),1,60);
    ExtraDataBeforeStarts = repmat(cvs(find(Vsec_after_gaps == 1)),1,60);

    ExtraDataAfterStops = reshape(ExtraDataAfterStops',length(igaps)*60,1);
    ExtraDataBeforeStarts = reshape(ExtraDataBeforeStarts',length(igaps)*60,1);

    % may bomb on these big matrices
    mat = [[Vsec; ExtraTimesAfterStops; ExtraTimesBeforeStarts] [cvs; ExtraDataAfterStops; ExtraDataBeforeStarts] ...
            [zeros(length(Vsec),1); ones(length([ExtraTimesAfterStops; ExtraTimesBeforeStarts]),1)]];
    mat = sortrows(mat);
    
    % alternative way to put together & sort data (rather than creating multi-column matrices):
    %Vsec = [Vsec; ExtraTimesAfterStops; ExtraTimesBeforeStarts];
    %cvs = [cvs; ExtraDataAfterStops; ExtraDataBeforeStarts];
    %[Vsec,i_sort] = sort(Vsec);
    %cvs = cvs(i_sort);

else  
    mat = [Vsec cvs];
    
end % if ~isempty...


% Smooth cumulative position data (2 sec hamming window) in deg/frame (60 frames/sec)
% (only roughly 2 sec - more accurate just to say "120-frame" hamming win. - esp. when
% there are timegaps in the data)
cvs_tsd=tsd(mat(:,1)*10000,mat(:,2));
%cvs_tsd=tsd(Vsec*10000,cvs); clear Vsec cvs;
Xs = SmoothTSD(cvs_tsd,120); clear cvs_tsd;

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

% calculate angular velocity
diff_sec=diff(Xs(:,1)/10000);
diff_data=diff(Xs(:,2));
dps=[0; diff_data./diff_sec];

dps_tsd = tsd(Xs(:,1),dps);