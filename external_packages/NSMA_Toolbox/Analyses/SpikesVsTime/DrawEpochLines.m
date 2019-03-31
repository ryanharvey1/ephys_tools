function DrawEpochLines(ymax,ymin,epochs)

% DrawEpochLines  Draws vertical lines representing epoch boundaries on graphs over time
%
% DrawEpochLines(ymax,ymin,epochs)
%
% INPUTS:  
%   epochs - a struct object with fields:
%       epochs.names = cell array of strings of epoch names 
%       epochs.intervals = cell array of 1x2 arrays with [start_ts  end_ts] (start
%           and end timestamps of each epoch) -- elements in cell array correspond
%           to epochs, listed in same sequence as in epochs.names
%   ymax - maximum y-axis value for lines to extend up to
%   ymin - minimum y-axis value for lines to extend down to
%   
% OUTPUTS:
%   (none)
%
% PL 2002, last modified '02 by PL


for iep = 1:length(epochs.names)
    t0 = epochs.intervals{iep}(1)/(60*10000);
    t1 = epochs.intervals{iep}(2)/(60*10000);
    tm = (t0+t1)/2;
    hhl = line([t0 t0], [ymax ymin]);
    set(hhl,'Color','r');
    hhl = line([t1 t1], [ymax ymin]);
    set(hhl,'Color','b');
    text(tm,ymax+(ymax-ymin)/40,epochs.names{iep});
end