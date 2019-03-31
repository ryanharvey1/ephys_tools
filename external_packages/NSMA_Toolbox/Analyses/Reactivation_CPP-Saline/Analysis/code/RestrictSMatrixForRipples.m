function SR = RestrictSMatrixForRipples(S,t0,t1)
% SR = RestrictTheSMatrix(S,t0,t1)
%
% Restricts the Spike Matrix S to the ripple intervals [t0, t1] (in timestamp 0.1 msec units),
% where t0 and t1 are arrays of timestamps of start and end times of eeg ripple activity.
% If t0 >= t1, an error is produced.
% 
% RS July 2004

[t0_m, t0_n] = size(t0);
[t1_m, t1_n] = size(t1);

if t0_m ~= t1_m
    %error('input arguments t0 and t1 must be the same size!');
    if (t0_m > t1_m) & (t0(end) > t1(end))              % if the start times list is longer than the end times list, and 
        t0(end) = [];                                   % the last start time is greater than the last end time, remove it
        
    elseif (t0_m < t1_m) & (t0(1) > t1(1))             % if the start times list is shorter than the end times list, and
        t1(1) = [];                                     % the first start time is greater than the first end time, remove it
            
    end 
end    

for i = 1:length(S)
    for j = 1:size(t0)
        if t0(j) >= t1(j)
            error('input argument t0 is requred to be smaller than t1!');
        else
            SR{i} = Restrict(S{i},t0(j),t1(j));
        end
    end
end