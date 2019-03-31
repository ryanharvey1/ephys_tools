function O = regularIntervals(start, stop, i_length, step, varargin);
% O = regularIntervals(start, stop, length, step) returns an intervalSet
% with intervals of constant length 
%
% INPUTS: 
% start, stop: of the period spanned by the intervalSet (intervalSet
% guaranteed to cover the period)
% length: length of each interval
% step: increment for each interval, if it's not present, defaults to
% length
% OUTPUT:
% O: the intervalSet
% OPTION:
% 'TimeUnits': the time unit in which the units are expressed (defaults to
% ts);
  
  
  
error(nargchk(3, 1000, nargin));
  
  if nargin == 3 
    step = i_length;
  end
  
  
  defined_options = dictArray( ...
      { { 'TimeUnits', {time_units('ts'), {'char', 'units'} } } });

  opt_varargin = varargin;

  getOpt;

  TimeUnits = time_units(TimeUnits);

  cnvrt = convert(TimeUnits, time_units('ts'));

  if cnvrt ~= 1
    start = start * cnvrt;
    stop = stop * cnvrt;
    i_length = i_length * cnvrt;
    step = step * cnvrt;
  end

  S = start:step:stop;

  if S(end) < stop - i_length
    S(end+1) = S(end) + step;
  end
  
  E = S+i_length;
  
  O = intervalSet(S, E);
  
  