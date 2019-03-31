function [t, f, S, clims] = plot_lfp_spectrogram(self, windowsize, windowinc, bandwidth, f_range, clims)
% Plots a spectrogram in jet colormap for current LFP
%
% Uses the Chronux toolbox.
%
% root.plot_lfp_spectrogram(lfp_ind, windowsize, windowsin, bandwidth,
% f_range);

    import CMBHOME.Utils.*
    
    AddChronuxPackage;
    
    if ~exist('windowsize', 'var')
        windowsize = 5;
    end
    
    if ~exist('windowinc', 'var')
        windowinc = 1;
    end
    
    if ~exist('f_range', 'var')
        f_range = [0 120];
    end
    
    if ~exist('bandwidth', 'var')
        bandwidth = 3;
    end
        
    signal = self.lfp.signal;
    
    if iscell(signal), signal = vertcat(signal{:}); end

    params.tapers = [bandwidth windowsize 1];
    params.fpass = f_range;
    params.Fs = self.lfp.fs;
    params.pad = 2;
    movingwin = [windowsize, windowinc];

    [S,t, f]=mtspecgramc(signal,movingwin,params);

    S = S*self.lfp.fs;
    
    if ~exist('clims', 'var'), clims = [0 max(max(S))/4]; end

    %imagesc(t, f, 10*log10(S)'), colorbar
    imagesc(t, f, S', clims), colorbar

end
