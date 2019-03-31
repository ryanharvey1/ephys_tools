function ThetaPowerVsPosition(V, thetapow_tsd, bincenters, horseshoe_flag)

% ThetaPowerVsPosition  Plots theta power and angular velocity vs. position (degrees) on circular track
%
% ThetaPowerVsPosition(V, thetapow_tsd, bincenters, horseshoe_flag)
%
% INPUTS:  
%   V - position tsd; angle on circular track in degrees [0..360) 
%   thetapow_tsd - tsd of timestamps and power at theta peaks, obtained from theta-filtered eeg
%   bincenters - centers of position bins in which to calculate avg. theta power and avg. velocity
%   horseshoe_flag - 1=animal runs in "horseshoe" pattern on circular track, where it runs once
%       around circle in one direction, turns and runs once around in reverse direction, etc.
%       0=animal runs regular laps around circle
%
% MRN 6/03


% THETA POWER VS. POSITION

%pow = Data(eeg_theta).^2;
%ts_eeg  = Range(eeg_theta,'ts');
%[ts_env,pow_env, PeaksIdx] = PowerEnvelope(ts_eeg,pow);
ts_env = range(thetapow_tsd,'ts');
pow_env = data(thetapow_tsd);

V_at_powts = data(restrict(V,ts_env));

% get mean theta power in each position bin
count = 1;
for i = bincenters;
    ix = find(V_at_powts >= i-2 & V_at_powts < i+2);
    if ~isempty(ix)
        pow_in_bins(count) = mean(pow_env(ix));
    else
        pow_in_bins(count) = NaN;
    end % if ~isempty...
    count = count+1;
end % for i

ix = find(pow_in_bins); 
imissing = find(isnan(pow_in_bins) == 1);
inotmissing = find(isnan(pow_in_bins) == 0);
pow_in_bins(imissing) = interp1(inotmissing,pow_in_bins(inotmissing),imissing);

%plot theta power vs. position
plot(bincenters,pow_in_bins,'k');  hold on;


% VELOCITY VS. POSITION
[dps_tsd] = GetSmoothedAngVelocity(V, horseshoe_flag);
dps = data(dps_tsd);

% get mean velocity in each position bin
Vdata = data(V);
count = 1;
for i = bincenters;
    ix = find(Vdata >= i-2 & Vdata < i+2);
    vel_in_bins(count) = mean(dps(ix));
    count = count+1;
end % for i

%plot velocity vs. position
vel_in_bins_rescaled = vel_in_bins * max(pow_in_bins)/max(vel_in_bins);
%vel_in_bins_rescaled = vel_in_bins;
plot(bincenters,vel_in_bins_rescaled,'g');

% axis labels for theta power & vel. plot
%xlabel('Position (degrees)');
%str(1) = {'\theta power &'};
%str(2) = {'velocity'};
%text(-.07,0,str,'Units','normalized','Rotation',90);
text(-.05,0,'\theta power','Color','k','Units','normalized','Rotation',90);
text(-.02,0,'Velocity','Color','g','Units','normalized','Rotation',90);
%ylabel(str);
if horseshoe_flag == 1
    axis([-360 360 0 max(vel_in_bins_rescaled)]);
else
    axis([0 360 0 max(vel_in_bins_rescaled)]);
end % if horseshoe
set(gca, 'XTickLabel', [], 'YTick', [], 'YTickLabel', []);