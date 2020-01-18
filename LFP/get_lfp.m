function data = get_lfp(lfpfile,data,varargin)
% get_lfp: loads Nlx csc lfp data, resamples it, and gets a few extra theta features  
%
% resample applies an antialiasing FIR lowpass filter to the signal and 
% compensates for the delay introduced by the filter
%
%   Input:
%           lfpfile: cell array containing file names of lfp files
%           data: ephys_tool's data structure
%           varargin:
%               Fold: lfp sample rate (32000)
%               Fnew: resampled lfp sample rate (1000)
%   Output:
%           data: 
%               data.lfp.ts
%               data.lfp.signal
%               data.lfp.theta
%               data.lfp.theta_phase
%               data.lfp.theta_amp
%               data.lfp.lfpsamplerate
%
% ryan h 2020

p = inputParser;
p.addParameter('Fold',32000);
p.addParameter('Fnew',1000);
p.parse(varargin{:});
Fold = p.Results.Fold;
Fnew = p.Results.Fnew;

% preallocate & get time stamps
TTnum=max(str2double(extractBetween(lfpfile,'CSC','.ncs')));
[ts] = Nlx2MatCSC(lfpfile{1},[1 0 0 0 0], 0, 1, [] );
signal = zeros(TTnum,(length(ts)*512)*Fnew/Fold);
theta=signal;
theta_phase=signal;
theta_amp=signal;

% resample time stamps & convert to sec for fma below
ts = interp1(linspace(1,length(signal),length(ts)), ts, 1:length(signal));
ts_sec=ts/10^6;
ts_sec=ts_sec-(ts_sec(1));
        
%account for non-integer fs
% get new fs exact
[fold, fnew] = rat(Fold./Fnew); 
Fnew = Fold.*(fnew./fold); 

% loop though each channel and resample lfp
fprintf('channel...');
for ii=1:length(lfpfile)
    fprintf(' %d ', ii);
    
    [~, filename] = fileparts(lfpfile{ii});
    ch = str2double(extractAfter(filename,'CSC'));

    [Samples]= Nlx2MatCSC(lfpfile{ii}, [0 0 0 0 1], 0, 1);

    signal(ch,:) = resample(Samples(:), fnew, fold);
        
    % filter for theta
    % normalized by the nyquist frequency
    Wn_theta=[4/(Fnew/2) 12/(Fnew/2)]; 
    [btheta,atheta]=butter(3,Wn_theta);
    theta(ch,:)=filtfilt(btheta,atheta,signal(ch,:));
    
    % FMA for phase & amp
    [phase,amplitude,~]=Phase([ts_sec',theta(ch,:)']);
    theta_phase(ch,:)=phase(:,2)';
    theta_amp(ch,:)=amplitude(:,2)';
end
fprintf('lfp loaded\n');

data.lfp.ts=ts;
data.lfp.signal=signal;
data.lfp.theta=theta;
data.lfp.theta_phase=theta_phase;
data.lfp.theta_amp=theta_amp;
data.lfp.lfpsamplerate=Fnew;
end