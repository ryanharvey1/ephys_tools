function O = set(O,varargin)

% Set intervalSet properties and return the updated object (TODO)
%  
%  TODO: what this function does???
%
% if 'start' is one property, 'stop' must be one, and viceversa

% batta 2001
% starting version

% Luke Sjulson, 2017-08: disabled this because I'm not entirely sure what
% it does or why it's necessary

warning('The set method is disabled');

% 
% 
% t = [];
% property_argin = varargin;
% set_start = 0;
% set_stop = 0;
% while length(property_argin) >= 2,
%     prop = property_argin{1};
%     val = property_argin{2};
%     property_argin = property_argin(3:end);
%     switch prop
%      case 'start'
%       set_start = 1;
%       start = val;
%      case 'stop'
%       set_stop = 1;
%       stop = val;
%      % here go more properties
%     end
% end
% 
% if xor(set_start, set_stop)
%   error('start and stop can only be set together');
% end
% 
% 
% % check start and stop, if necessary reshape and sort them, and update object
% if set_start
%   if ~isa(start, 'numeric')
%     error('start must be a numeric vector');
%   end
%   if all(size(start) ~= 1)
%     error('start must be a column or a row vector');
%   end
%   start = reshape(start, prod(size(start)), 1);
%   if ~isa(stop, 'numeric')
%     error('stop must be a numeric vector');
%   end
%   if all(size(stop) ~= 1)
%     error('stop must be a column or a row vector');
%   end
%   stop = reshape(stop, prod(size(stop)), 1);
% 
%   if length(start) ~= length(stop)
%     error('start and stop must be of the same length');
%   end
%   
%    if any(diff(start) < 0)
%     [start, ix] = sort(start);
%     stop = stop(ix);
%   end
%   
%   if any(start > stop)
%     error('Negative intervals are not allowed');
%   end
%   
%   O.start = start;
%   O.stop = stop;
% 
% end

