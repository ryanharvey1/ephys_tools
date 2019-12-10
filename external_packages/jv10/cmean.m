function y = cmean (x)
% CMEAN
% Circular equivalent to mean
%
% Ref: Statistical Analysis of Circular Data, N I Fisher
%
% --> www.paulbays.com

if any(abs(x)>pi), error('Input values must be in radians, range -PI to PI'); return; end

y = atan2(sum(sin(x)),sum(cos(x)));
