function [x,y] = FixPos(x,y,ts,max_allowed_flips)
% Fixes large jumps in position and interpolates missing values
%
% Takes all 0,0s and large jumps in position (greater than
% jitter_threshold=15 pixels/sample) that persist for less than
% max_allowed_flips (default = 5) and linearly interpolates the missing
% data. Smooths conservatively afterward, as well (convolution with a gaussian, standard
% deviation = 2 samples).
%

warning('off', 'MATLAB:interp1:NaNinY');

jitter_threshold=22;

bads = (x==0 | y==0);

x(bads) = NaN;
y(bads) = NaN;
ts_=ts;

if ~exist('max_allowed_flips', 'var')
    max_allowed_flips = 5; % samples
end

% [start,ends,~]=findgroups(isnan(x));
% max_allowed_flips=max(ends-start)+1;

flips = findOnsetsAndOffsets(isnan(x));

flips(:,2) = flips(:,2)+1;

flips = cat(1, 1, flips(:), find([0; sqrt(diff(x).^2 + diff(y).^2)]>jitter_threshold), length(x));

flips = sort(unique(flips)); % indeces of NaNs or jumps

flips = [flips(1:end-1), flips(2:end)];  % epochs formation

flips(:,2) = flips(:,2)-1; % adjust for diff shift

flips(flips(:,2)-flips(:,1)>max_allowed_flips-1,:) = [];

flips = mat2cell(flips, ones(size(flips,1),1),2); % convert to pairs corresponding to steps

flips = cellfun(@(c) c(1):c(2), flips, 'unif', 0); % convert to indices

x([flips{:}]) = []; % remove samples in ts and x
y([flips{:}]) = []; % remove samples in ts and x
ts([flips{:}]) = [];

x = interp1(ts, x, ts_);
y = interp1(ts, y, ts_);

x = ndnanfilter(x, normpdf(-6:6, 0, 2)', [], 1, {}, {}, 1); % conv with gaussian and ignore NaNs
y = ndnanfilter(y, normpdf(-3:3, 0, 2)', [], 1, {}, {}, 1);

% fill in remaining pesky NaN values
% XY=fillmissing([x,y],'pchip',1,'EndValues','nearest');
% x=XY(:,1);
% y=XY(:,2);
end

function [OnOffs] = findOnsetsAndOffsets(boolVec)
% Returns list of aligned start and stops of chunks of 1's in a vector of
% 1's and 0's.
%
% INPUTS
%  boolVec - vector of 1's and 0's
%
% OUTPUTS
%  startEnds - Nx2 list of indices of the first and last 1's for the N
%              contiguous blocks of 1's.
%
% % function [OnOffs] = findOnsetsAndOffsets(boolVec)


boolVec = boolVec(:)';

starts = find(diff(boolVec)==1);
ends = find(diff(boolVec)==-1);

% if starts out going correct speed, add 1 to starts
if boolVec(1)
    starts = [0 starts];
end

% if finishes going correct speed, add final value to ends
if boolVec(end)
    ends = [ends length(boolVec)];
end

OnOffs = [starts(:)+1 ends(:)];
end