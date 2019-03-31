function [ix_laps, nLaps, V_reparam] = LapCount(V, horseshoe_flag);

% LapCount  Breaks animal's path into laps on circle
%
% [ix_laps, nLaps, V_reparam] = LapCount(V, horseshoe_flag);
% 
% INPUTS:
%   V - position data tsd for circular track (degrees)
%   horseshoe_flag - 1=animal runs in "horseshoe" pattern where it runs a lap in one direction,
%       reverses direction in the next lap, reverses again in the next, etc.  
%       0(or no 2nd input)=animal runs regular laps around circle.
% OUTPUTS:
%   ix_laps - indices in reparameterized tsd where each lap ends (farthest point before animal turns 
%               around and begins next lap)
%   nLaps - number of laps run
%   V_reparam - (if horseshoe_flag=1) position data tsd reparameterized so that each lap is 
%       -360 -> 0 -> 360 degrees (a "lap" in this pattern is considered to be a full circle in one
%       direction followed by a full circle in the other)
%
% PL '01 & MN 6/03
%
% Recent Changes Made By Drew Maurer:
% There were (are?) problems with the horseshoe_flag code:
%
% if horseshoe_flag == 0... I don't know what Miriam was trying to acheive
% here, but it wasn't working. See Changes in code at lines 44-45 and 56-57

Vdata = data(V);
Vts = Range(V,'ts');
% take out dt zeros
dt=diff(Vts/10000);
Vts(find(dt == 0)+1,:)=[]; % remove 2nd frame of diff pair
Vdata(find(dt == 0)+1,:)=[];
V = tsd(Vts,Vdata);

if nargin < 2
    horseshoe_flag = 0;
end


if horseshoe_flag == 0  % rat runs regular laps around circle
    
    dV = diff(Vdata);
    ix_crossings = find((dV)< -300);
    %    ix_crossings = find(abs(dV)>300);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure; plot(dV)
%     hold on; plot(ix_crossings ,dV(ix_crossings), 'r.')
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    bad_dV_vals = dV(ix_crossings);
    corrected_dV_vals = -sign(bad_dV_vals).*(360-abs(bad_dV_vals));
    dV(ix_crossings) = corrected_dV_vals;
    cumV = cumsum(dV);
    
    %ix_laps = [ix_crossings(find(diff(cumV(ix_crossings)) > 300)); ix_crossings(end)];
    ix_laps = (find(diff(cumV) > 300));
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     figure; plot(cumV)
    %     hold on; plot(ix_laps, cumV(ix_laps), 'r.')
    %     figure; plot(dV);
    %     hold on; plot(ix_laps, dV(ix_laps), 'r.')
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(ix_crossings)
        if cumV(ix_crossings(1)) < 60    % (if animal briefly goes backwards and crosses 0-360 border
            ix_laps = ix_laps(2:end);    %  in beginning -- then "lap 1" is not a real lap - get rid of it)
        end % if cumV
    end

    nLaps = length(ix_laps);
    
    V_reparam = V;
    
end % if horseshoe_flag == 0


if horseshoe_flag == 1  % rat runs in "horseshoe" pattern on circular track
    
    % find degree measure of bottom of horseshoe where rat never runs
    bins = 1:359;
    count = hist(Vdata,bins);
    bottom = median(bins(find(count == min(count))));

    % rotate data
    crossing = mod(bottom+180,360);
    Vdata_rot = mod(Vdata-crossing,360);
    V_rot = tsd(Vts,Vdata_rot);

    % break path into laps
    dV = abs(diff(Vdata_rot));
    ix_crossings = find(dV>300);
    bad_dV_vals = dV(ix_crossings);
    corrected_dV_vals = 360-bad_dV_vals;
    dV(ix_crossings) = corrected_dV_vals;
    cumV = cumsum(dV);
    ix_laps = [ix_crossings(find(diff(cumV(ix_crossings)) > 60)); ix_crossings(end)];
    nLaps = length(ix_laps);

    % get real ix_laps (= last V index of each lap, = furthest point animal goes before turning around)
    nLminus1 = nLaps-1;
    ixlaps_temp = [];
    countlaps = 0;
    for iL = 1:nLminus1
        ts1 = Vts(ix_laps(iL));
        ts2 = Vts(ix_laps(iL+1));
        localV = restrict(V_rot,ts1,ts2);
        localVdata = data(localV);  localVrange = Range(localV,'ts');
        % check to make sure it's a real lap - that rat really gets to turnaround point in lap
        % if not, eliminate lap
        icheck = find(localVdata > 135 & localVdata < 225);
        if ~isempty(icheck)
            countlaps = countlaps+1;
            distance_from_bottom = abs(180-localVdata);
            i_bottom = find(distance_from_bottom == min(distance_from_bottom));
            ts_at_lap_end = localVrange(i_bottom(end));
            ixlaps_temp(countlaps) = find(Vts == ts_at_lap_end);
        else 
            nLaps = nLaps-1;
        end %if
    end
    ix_laps = ixlaps_temp;
    ix_laps(nLaps) = length(Vts);


    % change laps so that 1 lap = odd lap (-360-0) + even lap (0-360) -- 
    % 0 is bottom

    % get signed cumV
    dV_signed = diff(Vdata_rot);
    bad_dV_vals_signed = dV_signed(ix_crossings);
    corrected_dV_vals_signed = -sign(bad_dV_vals_signed).*(360-abs(bad_dV_vals_signed));
    dV_signed(ix_crossings) = corrected_dV_vals_signed;
    cumV_signed = [0; cumsum(dV_signed)];
    signlap1 = sign(cumV_signed(ix_laps(1)));

    ixLast = 0;
    new_ix_laps = zeros(floor(nLaps/2),1);
    for i=1:nLaps
        ix1=ixLast;
        ix2=ix_laps(i);  
        sign_this_lap = sign(cumV_signed(ix2)-cumV_signed(ix1+1));
        localcumv = cumV_signed(ix1+1:ix2);
        distthislap = localcumv-localcumv(1);  % start cum dist this lap from 0
        % make local cum dist increasing
        if sign_this_lap == -1 %if degrees decrease
            distthislap = -distthislap;
        end % if sign_this_lap == -1
        if sign_this_lap == signlap1 % if odd lap
            distfrombottom = abs(180-Vdata_rot(ix2));
            Vdata_rot(ix1+1:ix2) = distthislap-(distthislap(end)+distfrombottom);
        else % if even lap
            distfrombottom = abs(180-Vdata_rot(ix1+1));
            Vdata_rot(ix1+1:ix2) = distthislap+distfrombottom;
            new_ix_laps(i/2) = ix_laps(i);
        end %if sign_this_lap == signlap1
        ixLast = ix2;
    end %for i=...

    %ix_laps = new_ix_laps;
    if new_ix_laps(end) == length(Vts)
        ix_laps = new_ix_laps;
    else
        ix_laps = [new_ix_laps; length(Vts)];
    end %if
    nLaps = length(ix_laps);

    V_reparam = tsd(Vts, Vdata_rot);
    
end % if horsehoe_flag == 1