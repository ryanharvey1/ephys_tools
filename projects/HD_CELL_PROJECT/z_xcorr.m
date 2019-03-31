function [zhist] = z_xcorr(ts1, ts2)
% code to create z-scored cross-correlograms from Butler & Taube (2017)
% inputs: ts1 and ts2 (lists of spike times in microseconds)
% output: zhist (histogram with normalized counts in each timebin, 0.5 ms
% bins, -1000 to 1000 ms, 
% for zhist, 1st value is norm. count for -1000 ms lag, 2nd value is norm.
% count for -999.5 ms, ..., 2001st value is norm. count for 0 ms, etc.)

%% converts spike times to ms instead of microseconds
ts1 = (ts1/100);
ts2 = (ts2/100);

firstspiketime = min([ts1(:);ts2(:)]);
ts1adj = ts1-firstspiketime;
ts2adj = ts2-firstspiketime;
lastspiketime = max([ts1adj(:);ts2adj(:)]);

%converts spike time lists to histogram/counts of spike time differences
diffs = bsxfun(@minus,ts1adj,ts2adj');
small_diffs = diffs(abs(diffs) < 1000); % only keep spike differences up to 1000 ms
real_hist = hist(small_diffs,4001); % histogram with 0.5 ms bins

%% make distribution of 100 jittered/shuffled xcorrs
rep_num=0;
all_shuff_hist = nan(100,4001);
while rep_num<100
    jitter1 = ts1adj + randi([-10 10]); % jitter each spike train by +-10 ms
    jitter2 = ts2adj + randi([-10 10]);
    jitter_diffs = bsxfun(@minus,jitter1,jitter2');
    small_jitter_diffs = jitter_diffs(abs(jitter_diffs) < 1000); % only keep spike differences up to 1000 ms
    shuff_hist = hist(small_jitter_diffs,4001); % histogram with 0.5 ms bins
    all_shuff_hist(rep_num+1,:) = shuff_hist;
    rep_num=rep_num+1;
end

%compute average and sd for each shuffle lag bin to z-score actual histogram
mean_shuff = mean(all_shuff_hist);
sd_shuff = std(all_shuff_hist);
zhist = (real_hist-mean_shuff)./sd_shuff;