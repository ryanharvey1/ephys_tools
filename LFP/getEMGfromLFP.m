function  emg = getEMGfromLFP(lfp, varargin)

% based on Schomburg et al., 2014 and bz_EMGFromLFP. main differences are
% (1) independence from session info and xml file and (2) uses filtfilt
% instead of custom-designed. filteres lfp in the 300-600 Hz (changed to 300-499) band. bins
% the data according to the requested EMG sampling rate over sliding
% windows. window duration is determined by the ratio LFP / EMG sampling
% rate. Pearson's correlation coefficients is calculated at each bin for
% each pair of channels. the coefficients are averaged (summed and divided
% by the number of pairs). more information is availble in
% PreProcessing.docx
% 
% INPUT
%   lfp         matrix of samples (rows) x channels (columns)
%   fs          LFP sampling frequency {1250}
%   emgFs       EMG sampling frequency {40}
%   basepath    recording session path {pwd}
%   graphics    plot figure {1}.
%   saveVar     save variable {1}.
% 
% OUTPUT
%   emg         struct with fields:
%       data        EMG data (mean pairwise correlations)
%       timestamps  center of bin (window)
%       fs          EMG sampling frequency
% 
% TO DO LIST:
%   if emgFs < 1 then problem with graphics
% 
% 29 apr 19 LH. https://github.com/leoreh [slutskycode/extracellular in vivo/io/getEMGfromLFP.m]
%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'fs', 1250, @isnumeric);
addOptional(p, 'emgFs', 40, @isnumeric);
addOptional(p, 'basepath', pwd);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'saveVar', false, @islogical);

parse(p, varargin{:})
fs = p.Results.fs;
emgFs = p.Results.emgFs;
basepath = p.Results.basepath;
graphics = p.Results.graphics;
saveVar = p.Results.saveVar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nchans = size(lfp, 2);
chunksize = 20; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% option 1 (original)
% filter design. originally used by Schomberg. replaced with Matlab
% filtfilt because it is 10-times faster.
% tic
% maxf = floor(max([625 fs / 2]));       % adjust upper stopband to nyquist
% fband = [275, 300, maxf - 25, maxf];
% f  = fdesign.bandpass(fband(1), fband(2), fband(3), fband(4),...
%     60, 1, 60, fs);
% filt = design(f, 'butter', 'MatchExactly', 'passband');
% 
% % filter in both the forward and reverse directions
% filt_sig = zeros(size(lfp));
% for i = 1 : size(lfp, 2)
%     filt_sig(:, i) = filter(filt, lfp(:, i));
%     filt_sig(:, i) = filter(filt, filt_sig(end:-1:1, i));
%     filt_sig(:, i) = filt_sig(end:-1:1, i);
% end
% toc

%%% option 2
filtered = filterLFP(lfp, 'fs', fs, 'type', 'butter', 'passband', [300 600],...
    'graphics', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% correlation coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

win_dur = 1 / emgFs;
win_samps = round(win_dur * fs);
win_idx = -win_samps : win_samps;
timestamps = (1 + win_idx(end)) : win_samps : (size(lfp, 1) - win_idx(end));
nbins = length(timestamps);

data = zeros(nbins, 1);
for j = 1 : nchans
    for k = (j + 1) : nchans
        
        msg = sprintf('calculating correlations between channels %d and %d', j, k);
        disp(msg);
        
        c1 = [];
        c2 = [];
        bidx = 0;
        bstart = 1;
        for i = timestamps
            bidx = bidx + 1;
            s1 = filtered(i + win_idx, j);
            s2 = filtered(i + win_idx, k);
            c1 = cat(2, c1, s1);
            c2 = cat(2, c2, s2);
            if size(c1, 2) == chunksize || i == timestamps(end)
                bend = bidx;
                tmp = corr(c1, c2);
                tmp = diag(tmp);
                data(bstart : bend) = data(bstart : bend) + tmp;
                c1 = [];
                c2 = [];
                bstart = bidx + 1;
            end
        end
    end
end

data = data / (nchans * (nchans - 1) / 2); % divide with number of pairs

emg.data = data;
emg.timestamps = timestamps' / fs;
emg.fs = emgFs; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plots maximum and minimum EMG activity 
if graphics
    figure
    interval = 1.5;
    
    subplot(1, 2, 1)
    [~, idx] = max(emg.data);
    lfpidx = idx / emg.fs * fs;
    idx = idx - interval * emg.fs : idx + interval * emg.fs;
    lfpidx = round(lfpidx - interval * fs : lfpidx + interval * fs); 
    yfactor = max(lfp(lfpidx, 1)) - min(lfp(lfpidx, 1));
        
    yyaxis left
    plot(idx / emg.fs / 60, emg.data(idx), 'LineWidth', 3)
    ylabel('EMG (mean pairwise correlation)')
    hold on
    yyaxis right
    for i = 1 : nchans
        plot(lfpidx / fs / 60, lfp(lfpidx, i) + i * yfactor, '-k')
    end
    ylabel('LFP raw', 'Color', 'k')
    set(gca,'YTick',[])
    xlabel('Time [m]')
    title('maximum correlation')
    axis tight
    
    subplot(1, 2, 2)
    [~, idx] = min(emg.data);
    lfpidx = idx / emg.fs * fs;
    idx = idx - interval * emg.fs : idx + interval * emg.fs;
    lfpidx = round(lfpidx - interval * fs : lfpidx + interval * fs); 
    yfactor = max(lfp(lfpidx, 1)) - min(lfp(lfpidx, 1));
        
    yyaxis left
    plot(idx / emg.fs / 60, emg.data(idx), 'LineWidth', 3)
    ylabel('EMG (mean pairwise correlation)')
    hold on
    yyaxis right
    for i = 1 : nchans
        plot(lfpidx / fs / 60, lfp(lfpidx, i) + i * yfactor, '-k')
    end
    ylabel('LFP raw', 'Color', 'k')
    set(gca,'YTick',[])
    xlabel('Time [m]')
    axis tight
    title('minimum correlation')
    sgtitle('EMG from LFP')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if saveVar   
    [~, filename] = fileparts(basepath);
    save([basepath, '\', filename, '.EMGfromLFP.mat'], 'emg')
end

end