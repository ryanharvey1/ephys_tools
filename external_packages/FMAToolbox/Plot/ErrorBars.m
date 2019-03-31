%ErrorBars - Plot mean±SEM.
%
% Plot mean±SEM for repeated observations.
%
%  USAGE
%
%    p = ErrorBars(X,<options>)
%
%    X              the matrix to plot; each column is an observation
%    <options>      options for function
%

% Copyright (C) 2004-2006 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function p = ErrorBars(X,varargin)

if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help ErrorBars">ErrorBars</a>'' for details).');
end

p = errorbar(mean(X,2),sem(X));
set(p,varargin{:});
