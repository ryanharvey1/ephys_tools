function y = cstd (x)
% CSTD
% Circular equivalent to standard deviation
%
% Ref: Statistical Analysis of Circular Data, N I Fisher
%
% --> www.paulbays.com

if any(abs(x)>pi), error('Input values must be in radians, range -PI to PI'); return; end

if size(x,1)==1, x = x'; end

R = sqrt(sum(sin(x)).^2+sum(cos(x)).^2) / size(x,1);
y = sqrt(-2 * log(R));