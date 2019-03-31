function RasterPlot(S, varargin)
% RasterPlot(S, height, bar, tstart, tend)
% 
% Plot spike trains as rasters.
%
% INPUT:
% S: cell array of ts objects, spike trains
%
% OPTIONAL INPUTS: 
% height: the height allocated to each cell (includes spacing), default 1
% bar: the fraction of height occupied by a cell, default 0.8
% tstart, tend: the beginning of end of time interval to be plot, default
% the whole thing
  
% batta 2001 under construction
  
height = 1;
if nargin >= 2
  height = varargin{1};
end

bar_fraction = 0.8;
if nargin >= 3
  bar_fracion = varargin{2};
end

tstart = 0;
if nargin >= 4 
  tstart = varargin{3};
end

tend = -1;
if nargin >= 5
  tend = varargin{4};
end

if tstart ~= 0 & tend < 0 % then you have to find the max in S
  for i = 1:length(S)
    mm(i) = EndTime(S{i})
  end
  tend = max(mm);
end

if tstart > 0
  for i = 1:length(S)
    Si{i} = Restrict(S{i}, tstart, tend);
  end
else 
  Si = S;
end


  
for i = 1:length(Si)
  sp = Data(Si{i});
  sx = [sp sp repmat(NaN, length(sp), 1)];
  sy = repmat([(i*height) (i*height + height *bar_fraction) NaN], length(sp), 1);
  sx = reshape(sx', 1, length(sp)*3);
  sy = reshape(sy', 1, length(sp)*3);
  
  line(sx, sy, 'Color', 'k');
  hold on
end

