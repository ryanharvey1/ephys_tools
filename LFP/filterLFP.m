function  filtered = filterLFP(raw, varargin)

% applies filtfilt using different parameters
% 
% INPUT
%   raw         matrix of samples (rows) x channels (columns). must be
%               double
%   fs          sampling frequency
%   passband    pass frequency range. [0 X] for low-pass, [X inf] for highpass
%   stopband    stop frequency range
%   order       filter order. if cheby2 order is number of cycles {4}
%               if fir1 order in samples
%   type        filter type {'cheby2'}, 'fir1' or 'butter'
%   rs          filter ripple {20}
%   dataOnly    logical. return only data or also struct
%   basepath    recording session path {pwd}
%   graphics    plot figure {1}.
%   saveVar     save variable {1}.
% 
% OUTPUT
%   if dataOnly
%       filtered    vector with filtered data
%   else
%       filtered    struct with fields:
%           data        filtered data
%           amp         amplitude, calculated via hilbert transform
%           phase       phase angle, calculated via hilbert transform
%           <params>    as in input
% 
% 24 apr 19 LH. https://github.com/leoreh 
%       29 apr 19       added dataOnly

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'passband', [], @isnumeric);
addOptional(p, 'stopband', [], @isnumeric);
addOptional(p, 'fs', 1250, @isnumeric);
addOptional(p, 'order', 4, @isnumeric);
addOptional(p, 'rs', 20, @isnumeric);
addOptional(p, 'type', 'cheby2', @ischar);
addOptional(p, 'dataOnly', true, @islogical);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'saveVar', false, @islogical);

parse(p, varargin{:})
basepath = p.Results.basepath;
fs = p.Results.fs;
passband = p.Results.passband;
stopband = p.Results.stopband;
order = p.Results.order;
rs = p.Results.rs;
type = p.Results.type;
dataOnly = p.Results.dataOnly;
graphics = p.Results.graphics;
saveVar = p.Results.saveVar;

% params
nyquist = fs / 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch(type)
	case 'cheby2'
		if ~isempty(passband)
			if passband(1) == 0
				[b, a] = cheby2(order, rs, passband(2) / nyquist, 'low');
            elseif passband(2) == inf
                [b, a] = cheby2(order, rs, passband(1) / nyquist, 'high');
			else
				[b, a] = cheby2(order, rs, passband / nyquist);
			end
		else
			[b a] = cheby2(order, rs, stopband / nyquist, 'stop');
        end
	case 'fir1'
		if ~isempty(passband)
			if passband(1) == 0
                filt_order = round(order * 2 * nyquist ./ passband(2));    
				[b a] = fir1(filt_order, passband(2) / nyquist, 'low');
            elseif passband(2) == inf
                filt_order = round(order * 2 * nyquist ./ passband(1));    
				[b a] = fir1(filt_order, passband(1) / nyquist, 'high');
            else
                filt_order = round(order * 2 * nyquist ./ passband(1));  
				[b a] = fir1(filt_order, passband / nyquist);
			end
        else
            filt_order = round(order * 2 * nyquist ./ stopband(1));
			[b a] = fir1(filt_order, stopband / nyquist, 'stop');
        end
    case 'butter'
        if ~isempty(passband)
            [b a] = butter(order, [passband(1) / nyquist passband(2) / nyquist], 'bandpass');
        else
            [b a] = butter(order, stopband(1) / nyquist, 'stop');
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : size(raw, 2)
    data(:, i) = filtfilt(b, a, raw(:, i));
end

if ~dataOnly
    filtered.data = data;
    for i = 1 : size(raw, 2)
        hilb = hilbert(filtered.data(:, i));
        filtered.amp(:, i) = abs(hilb);
        filtered.phase(:, i) = angle(hilb);
    end
    
    filtered.type = type;
    filtered.passband = passband;
    filtered.stopband = stopband;
    filtered.order = order;
    filtered.fs = fs;
else
    filtered = data;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    t = (1 : length(data)) / fs;
    h = figure;
    p1 = plot(t, data(:, 1), 'LineWidth', 2);
    hold on
    p2 = plot(t, raw(:, 1));
    p2.Color(4) = 0.8;
    legend('Filtered', 'Raw')
    axis tight
    xlabel('Time [s]')
    xlim([1 2])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveVar   
    [~, filename] = fileparts(basepath);
    save([basepath, '\', filename, '.filtered.mat'], 'filtered')
end

end

% EOF