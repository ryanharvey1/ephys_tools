%function [y b a]=ButFilter(x, n, wn, flag)
% x - signal, n - order, wn - freq.range of filter, flag - 'bandpass',
% 'low', or 'high' in normalized units: Hz/(SamplingRate/2)
function [y b a] =ButFilter(x, varargin)
[n, wn, flag] = DefaultArgs(varargin, {20, 50/625,'low'});
[b,a] = butter(n,wn,flag);
nCh = size(x,2);
%if size(x,1) == nCh
%    x=x';
%end

for i=1:nCh
    y(:,i)=filtfilt(b,a, x(:,i));
end

