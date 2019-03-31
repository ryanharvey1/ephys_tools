function [t, p, ci, ndof, h] = SlopeDiff_TTest(s1,sig1,n1,s2,sig2,n2, alpha, tail)
%
%  [t, p, ci, ndof, h] = SlopeDiff_TTest(s1,sig1,n1,s2,sig2,n2, alpha, tail)
%
%  given 2 regression slopes s1, s2 with two slope standard-errors sig1, sig2 and the 2 sample sizes n1,n2 these regression slopes are
%  based on, compute the two-tailed (for tail = 0, default) t-statistic and significance level p of the hypothesis that the difference of slopes s1-s2 = 0 (i.e the null-hypothesis h=0
%  that the slopes are the same).
%
% INPUT: 
%    s1, s2 ... regression slopes as returned by regress
%    sig1, sig2 ... STD of the sampling distribtion of s1 and s2 (see below)
%    n1, n2 ...  sample sizes (number of points) on which the regressions slopes are based on ( e.g n1 = length(x) in the regression model below)
%    alpha  .... confidence level for accepting the null-hypotheses (the default alpha = 0.05 means 95% confidence level for REJECTING it) 
%    tail   ... 0 = 2-tailed t-test,  (if h=1 "s1 and s2 are (1-alpha)-significantly different") 
%              -1 = left tailed test (if h=1 "s1 is (1-alpha)-significantly less than s2")
%              +1 = right tailed test (if h=1 " s1 is (1-alpha)-significantly greater than s2")
%
%
% OUTPUT:
%    t  ... value of the t-statistic t = s1-s2/s_D, where s_D is the 'pooled std of the mean' s_D 
%    p  ... p-value (significance; prob. that an experiment gives that t-value by chance)
%    ci ... vector [c1, c2] of the 1-alpha confidence interval of the slope difference s1-s2
%    ndof ... Number of degree of freedoms (ndof= n1+ n2 -2)
%    h  ... accept (h=0) or reject (h=1) the null-hypotheses that the slopes s1 and s2 are equal.
%  
% Since the regression slopes of each indiviual regression line have a gaussian sampling distribution 
% (cf. 'Applied Linear Regression Models' by J. Neter, W. Wasserman and M.H. Kutner)
% so are their differences and we can do a simple 2-tail t-test.
%
% The standard errors of the mean ( = std of the sampling distributions) sig1 and sig2 of s1 and s2 can be calculated as follows:
% a) approx: 
%    if the n1 and n2 >> 1, the t-distribution is almost gaussian and you can estimate the
%    sig1 and sig2 by e.g sig1= abs(c1-c2)/2, where c1 and c2 are the bounds of a 95.45% (=approx 2*sig1)
%    confidence interval as returned by the matlab regress function.
%
% b) exact: 
%    sig1^2 = sum(r.^2)/((n1-2)*sum((x-x_avg).^2))     (standard error estimator of regression slope; formula 3.9 in 
%           (and similar for sig2^2) where
%    r ... n1 vector of residuals (as returned by the regress function)
%    x ... n1 vector of x-values that enter into the linear regression model y = s1*x + s0 + r where r is normal(0,sigma) distributed
%          (s1 and s0 are the regression slopes and intercepts respectively, as returned by regress)
%    x_avg = mean(x)
% 
%
% needs the stats toolbox (for tinv and tcdf)
%
% Literature:
%  John Neter, William Wasserman and Michael H. Kutner 'Applied Linear Regression Models 2nd Ed', IRWIN (1989)
%  William Press et al. 'Numerical Recipies in C, 2nd Ed.' Cambridge (1992)
%
%  Peter Lipa  April 12, 2004


% check inputs and handle default arguments
if nargin < 7 | isempty(alpha)
    alpha = 0.05;
end

if (prod(size(alpha))>1), error('alpha must be a scalar.'); end
if (alpha<=0 | alpha>=1), error('alpha must be between 0 and 1'); end

if nargin < 8, 
    tail = 0; 
end 

if isempty(find([-1,0,1]-tail)), error('tail must be one of -1, 0, 1 exactly'); end
    

% numerator of t-statistic
xmean = s1 - s2;

% denominator of t-statistic: pooled variance s_D (formula 14.2.1 in Press et al. 'Numerical Recipies')
ss = (n1*sig1)^2 + (n2*sig2)^2;
s_D = sqrt((ss/(n1+n2-2))*(1/n1 + 1/n2));
%s_D = sqrt(sig1*2/n1 + sig2*2/n2);                % formula 14.2.3 in Press et al. for unequal variance t-test  (not implemented Yet)

% t-statistic
t = xmean / s_D;

% number of degree of freedoms
ndof = n1+n2-2;

% p-value from Student cumulative distribution function at t
p = tcdf(t,ndof);          % the p-value just found is for the  tail = -1 test
if (tail == 0)
    p = 2 * min(p, 1-p);       % p-value for the 2-tailed t-test
    % 1-alpha confidence interval for the slope difference
    crit = tinv(1-alpha/2,ndof) * s_D; 
    ci = [(xmean - crit) (xmean + crit)];
else
    crit = tinv(1-alpha,ndof) * s_D;   % for tail = +-1
    if tail == 1
        p = 1 - p;
        ci = [(xmean - crit), Inf];
    else
        ci = [-Inf, (xmean + crit)];
    end    
end

% Determine if the actual significance exceeds the desired significance
h = 0;
if p <= alpha, 
    h = 1; 
end 

if isnan(p), 
    h = NaN; 
end



