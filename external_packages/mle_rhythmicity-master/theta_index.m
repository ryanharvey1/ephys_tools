function [thetaindex,p,f,F,theta_index_ci,lowsig,highsig,blow,bhigh] = theta_index(x_,varargin)
%THETA_INDEX Calculates the theta index and confidence intervals
% Calculates the theta-index as in Yartsev et al, 2011, Nature. See also:
% Boccara et al, 2010, Nature Neurosci; Langston et al, 2010 Science; Wills
% et al, 2010, Science; and Deshmukh et al, 2010, J. Neurophysiol.
%
% Using a model of the autocorrelogram as a series of Poisson counts, also
% calculates confidence intervals on the theta index.
%
% INPUT
%   x: Data. Can either be spike timestamps, a cell array of lags, or a
%       histogram of lag counts.
% PARAMETERS
%   calcp (true): Bootstaps the null distribution by shuffling lags by 1/2
%       of a cycle of the peak frequency detected.  If false, does not run
%       and p=NaN, but runs much faster.
%   pmethod (bootstrap): Can be "bootstrap" or "shuffle". If "shuffle", x
%       must be spike timestamps and generates the null distribution by
%       jittering spike time as specified by shufflet.  If "bootstrap",
%       jitters spike lags within the AC window by the amount specified by
%       shufflet.
%   shufflet ('adapt'): If 'adapt', jitters spike/lag times by +/- (1/(f-1)/2)
%       seconds (see output).  Can also be a set number (i.e. 10 seconds).
%   t (0:0.01:0.5): The edges of the lag bins.
%   ishist (false): Set to true if input data is a histogram.
%   calcci (true): Calculates the region of confidence (0.95) around the
%       autocorrelogram for the true underlying rates, and finds the least
%       and most theta-rhythmic (by the theta index) rate profiles in this
%       range. If false, does not run and
%       theta_index_ci,lowsig,highsig,blow and bhigh are all NaN, but runs
%       much faster
% RETURNS
%   thetaindex: The theta index - the ratio of the average power within 1
%       Hz of the peak in the theta (5-11 Hz) range and the average power
%       in the whole spectrum of the zero-mean spike time autocorrelation
%   p: The boostrapped significance of the autocorrelation
%   f: The peak theta frequency detected
%   F: The bootstrapped null distribution
%   theta_index_ci: The confidence intervals on the theta index.
%   lowsig: The rate profile found with the lowest theta index.
%   highsig: The rate profile found with the highest theta index.
%   blow: The low side of the 95% confidence range for the rate profile.
%   bhigh: The high side of the 95% confidence range for the rate profile.
%
% Copyright (c) 2014, Trustees of Boston University
% All rights reserved.
%
% This file is part of mle_rhythmicity
%
% This code has been freely distributed by the authors under the BSD 
% licence (http://opensource.org/licenses/BSD-2-Clause). If used or
% modified, we would appreciate it if you cited our paper:
%
% Climer, J. R., DiTullio, R., Newman, E. L., Hasselmo, M. E., Eden, U. T. 
% (2014), Examination of rhythmicity of extracellularly recorded neurons in
% the entorhinal cortex. Hippocampus, Epub ahead of print. doi:
% 10.1002/hipo.22383.

%% Parse input
p = inputParser;
p.addParamValue('calcp',true);
p.addParamValue('ishist',false);
p.addParamValue('calcci',true);
p.addParamValue('t',0:0.01:0.5);
p.addParamValue('pmethod','bootstrap');
p.addParamValue('shufflet',10);
p.parse(varargin{:});
for j = fields(p.Results)'
    eval([j{1} ' = p.Results.' j{1} ';']);
end

T = max(t);


if ishist
    b = x_(:);
elseif iscell(x_)
    x = x_(:);
    x_ = [];
    try
        b = histc(cat(1,x{:}),t);
    catch err
        b = histc(cat(2,x{:}),t);
    end
    b = b(1:end-1);
else
    spk_ts = x_(:);
    x_ = x_(:)';
    x = arrayfun(@(x)spk_ts(spk_ts>x&spk_ts<=x+T)-x,spk_ts,'UniformOutput',false);
    b = histc(cat(1,x{:}),t);
    b = b(1:end-1);
end

b = b(:);

if numel(b)==numel(t)-1
    b = [flipud(b);max(b(:));b];
end


% Calculate theta index
Fs = 0.01^-1;
NFFT = 2^16;
f = Fs/2*linspace(0,1,NFFT/2+1);
Y = fft(b-mean(b),NFFT);
Y = abs(Y(1:NFFT/2+1)).^2;
Y = conv(Y,ones(1,round(2/mode(diff(f)))),'same')';
[~,pk] = nanmax((f(:)>5&f(:)<11).*Y(:));
pk0 = pk;

thetaindex = mean(Y(abs(f-f(pk))<=1))/mean(Y);
p = NaN;
F = NaN;
theta_index_ci = NaN;

%% Calculate confidence intervals on theta index
if calcci
    if (isempty(x_)||~iscell(x_))&&iscell(x)
        n = numel(x);
    elseif iscell(x_)
        n = numel(x_);
    end
    
    b2 = b(numel(t)+1:end);
    blow = b2*0;
    blow(b2~=0) = arrayfun(@(b2)fzero(@(mu)poisscdf(b2,mu*n)-1+0.05/2,[0 10*b2/n]),b2(b2~=0))*n;
    bhigh = arrayfun(@(b2)fzero(@(mu)poisscdf(b2,n*mu)-0.05,[0 10*max(b2,1)/n]),b2)*n;
    
    theta_index_ci = [0 0];
    
    % Find underlying rate profile with lowest theta index
    [lowsig,theta_index_ci(1)] = fmincon(...
        @(x)confun(x,bhigh,blow)...
        ,b2...
        ,[],[],[],[]...
        ,blow...
        ,bhigh...
        ,[]...
        ,optimset('Algorithm','active-set','TolFun',1e-4,'TolX',1e-4)...
        );
    
    % Find underlying rate profile with highest theta index
    [highsig,theta_index_ci(2)] = fmincon(...
        @(x)-confun(x,bhigh,blow)...
        ,b2...
        ,[],[],[],[]...
        ,blow...
        ,bhigh...
        ,[]...
        ,optimset('Algorithm','active-set','TolFun',1e-3,'TolX',0.01)...
        );
    
    theta_index_ci(2) = -theta_index_ci(2);
end

%% Bootstrap: uniformly jitter lag by 1/2 cycle and recalculate
N = 1000;

if calcp
    % Calc shufflet if nessesary
    if isequal(shufflet,'adapt')
        shufflet = 1/(f(pk0)-1)/2;
    end
    
    F = NaN(N,1);
    
    switch pmethod
        case 'bootstrap'
            for j=1:N
                % jitter
                try
                    shuffled_b = abs(unifrnd(cat(1,x{:})-shufflet,cat(1,x{:})+shufflet));
                catch err
                    shuffled_b = abs(unifrnd(cat(2,x{:})-shufflet,cat(2,x{:})+shufflet));
                end
                
                % Make jittered histogram
                shuffled_b(shuffled_b>T) = 2*T-shuffled_b(shuffled_b>T);
                shuffled_b(shuffled_b<=0) = -shuffled_b(shuffled_b<=0);
                shuffled_b = histc(shuffled_b,t);
                shuffled_b = shuffled_b(1:end-1);
                shuffled_b = shuffled_b(:);
                shuffled_b = [flipud(shuffled_b);max(shuffled_b);shuffled_b];
                
                % Calculate theta index
                Y = fft(shuffled_b-mean(shuffled_b),NFFT);
                Y = abs(Y(1:NFFT/2+1)).^2;
                Y = smooth(Y,2/mode(diff(f)));
                [~,pk] = nanmax((f>5&f<11).*Y');
                F(j) = mean(Y(abs(f-f(pk))<=1))/mean(Y);
            end
        case 'shuffle'
            
            for i=1:N
                % Calculate theta index
                F(i) = theta_index(spk_ts(:)+rand(numel(spk_ts),1)*shufflet*2-shufflet,'calcp',false,'calcci',false,'t',t);
            end
        otherwise
            warning('theta_index:BadPMETHOD',['pmethod ''' pmethod ''' not implemented.']);
    end
    p = sum(F>thetaindex)/N;
end

f = f(pk0);

end

function [ c ] = confun( x, high, low )
c = theta_index(x,'calcp',false,'ishist',true,'calcci',false);
end