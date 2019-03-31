function [p,U2_obs,U2_H0]=watsons_U2_perm_test(A1,A2,N)
%
% [p,U2_obs,U2_H0]=watsons_U2_perm_test(A1,A2,N)
%
% Computes the probability that two samples of circular data (angles) come
% from the same distribution, or from distributions with the same direction
% (null hypothesis), using a permutation test and Watson's U2 statistic.
% This function calls watsons_U2 to compute that statistic. The p values
% obtained with this code have been verified using the numerical examples
% 27.10 and 27.11 from Zar (1999).
%
% Inputs:
% A1, A2:   vectors containing angles (in degrees or radians, the unit does
%           not matter)
% N:        integer scalar containing the number of permutations to perform
%
% Outputs:
% p:        probability that the two samples come from the same
%           distribution
% U2_obs:   observed value of Watson's U2 statistic
% U2_H0:    distribution of Watson's U2 statistic under the null hypothesis
%
% References:
% Zar JH. Biostatistical Analysis. 4th ed., 1999.
% Chapter 27: Circular Distributions: Hypothesis Testing.
% Upper Saddle River, NJ: Prentice Hall.
%
% pierre.megevand@gmail.com
%
% See also WATSONS_U2

U2_obs=watsons_U2(A1,A2);
U2_H0=zeros(N,1);
try
    A=[A1;A2];
catch
    A=[A1 A2];
end
n1=numel(A1);
n2=numel(A2);
n12=n1+n2;

for i=1:N
    a=A(randperm(n12));
    a1=a(1:n1);
    a2=a(n1+1:end);
    U2_H0(i)=watsons_U2(a1,a2);
end

p=nanmean(U2_H0>=U2_obs);
