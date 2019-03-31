function h = fit(varargin)

% $Revision: 1.17.2.2 $  $Date: 2005/03/07 17:24:42 $
% Copyright 2000-2005 The MathWorks, Inc.


% We may be asked to create an empty object not connected to the database
toconnect = 1;
if nargin==1 && isequal(varargin{1},'disconnected')
   varargin(1) = [];
   toconnect = 0;
end

% The meat of the constructor is in the helper function, constructorhelper.
% See the comments in that file for more information.
h = constructorhelper(cftool.fit,varargin{:});

if toconnect
   connect(h,getfitdb,'up');
end
