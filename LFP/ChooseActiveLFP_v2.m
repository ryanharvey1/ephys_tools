function [f_peak, active_lfp_ind] = ChooseActiveLFP_v2(self)

% active_lfp_ind = ChooseActive(root)

        for k = 1:size(self.signal,1)
           [f_peak(k), pr(k)] = calc_spectrum(self.signal(k,:),self.lfpsamplerate);
        end
        
        [~, active_lfp_ind] = max(pr);
        f_peak=f_peak(active_lfp_ind);
        
end


function [f_peak,power_ratio] = calc_spectrum(lfp,fs)

    TBP = 10;
    n_tapers = 9;
    signal = lfp;
    
    f_low = 0; % Hz
    f_high = 20; % Hz
    f_range = [f_low, f_high];
    bandwidth = .25; %Hz
    windowsize = 10; % seconds
    windowinc = 2; % seconds
    TBP = 10; 
    tapers = 9;
    
    if iscell(signal), signal = vertcat(signal{:}); end
        
    params.tapers = [TBP, n_tapers];
    params.fpass = f_range;
    params.Fs = fs;
    params.pad = 2;

    [S,f]=mtspectrumc(signal,params);

    stdf = .2/mean(diff(f)); % elements in a .2 hz std
    
    kernel = pdf('Normal', -3*stdf:3*stdf, 0, stdf); % smooth with gaussian kernel
    
    S = conv(S, kernel, 'same'); % plot(f, S)
    
    [xmax,imax] = extrema(S); % find local maxima and see if any lie within theta

    ftmp = f(imax);
    stmp = S(imax);

    [a_peak,ind] = max(stmp(ftmp > 7 & ftmp < 11)); % intrinsic frequency

    ftmp = ftmp(ftmp > 7 & ftmp < 11);
    
    if ~isempty(ind)                                    % if theta peak was found
        f_peak = ftmp(ind)';     
        a_peak_av = mean(S(f>f_peak-1 & f<f_peak+1));
        peak = S(f==f_peak);
    else                                                % if theta peak was not found
        f_peak = 0;
        peak = 0;
        a_peak_av = 0;
    end

    power_ratio = a_peak_av / mean(S);
    
end
