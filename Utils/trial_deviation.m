function d = trial_deviation(trials)
% trial_deviation: computes standard deviation across spatial bins over
% trials
%
% Input
%   trials: n by m matrix where n = trials, m = spatial bins
%
% Output
%   d: mean standard deviation over trials across spatial bins
%
trials(isnan(trials)) = 0;

norm_t = trials - min(trials(:));
norm_t = norm_t ./ max(norm_t(:));

d = mean(nanstd(norm_t,0,1));
end

