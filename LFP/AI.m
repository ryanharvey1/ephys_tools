data=load('d:\Users\BClarkLab\Desktop\HPCatn02_S20180707144931');
signal=data.lfp.signal;

%Broad Bandpass filter for determine local minima & maxima
for ch=1:size(signal,1)
broad_signal(ch,:) = BandpassFilter(signal(ch,:), data.lfp.lfpsamplerate, [1, 80]); 
end

%Find correlated channels
corMat=corrcoef(broad_signal');

%Check Correlations - flag greater than .7
check=corMat >=.7;

%Compute median signal
medSig=median(broad_signal,1);

%Find local minima 
[troughIdx,~] =islocalmin(medSig,'MinSeparation',(data.lfp.lfpsamplerate/12));
[peakIdx,~] =islocalmax(medSig,'MinSeparation',(data.lfp.lfpsamplerate/12));


% 
% ascend=
% descend=fing(medSig(trough):medSig(peak));


figure;
plot(medSig(1:1000));hold on
plot(P(1:1000));hold on

plot()

%Extra phase via linear interpolation between local min/max


%################# LOCAL FUNCTIONS ########################################

    function signal_filtered = BandpassFilter(signal, Fs, Fpass)
        % Takes 'signal' and bandpasses it to Fpass frequencies
        %
        % Arguments
        %
        % signal - arbitrary signal to be filtered
        % Fs - sample frequency
        % Fpass - 1x2 vector of F_low and F_high indicating passband
        %
        % signal_filtered = BandpassFilter(signal, Fs, Fpass)
        
        Wn_theta = [Fpass(1)/(Fs/2) Fpass(2)/(Fs/2)]; % normalized by the nyquist frequency
        
        [btheta,atheta] = butter(3,Wn_theta);
        
        signal_filtered = filtfilt(btheta,atheta,signal);
    end