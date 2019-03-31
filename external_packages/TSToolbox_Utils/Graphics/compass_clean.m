function hc = compass(varargin)
%COMPASS Compass plot.
%   COMPASS(U,V) draws a graph that displays the vectors with
%   components (U,V) as arrows emanating from the origin.
%
%   COMPASS(Z) is equivalent to COMPASS(REAL(Z),IMAG(Z)). 
%
%   COMPASS(U,V,LINESPEC) and COMPASS(Z,LINESPEC) uses the line
%   specification LINESPEC (see PLOT for possibilities).
%
%   COMPASS(AX,...) plots into AX instead of GCA.
%
%   H = COMPASS(...) returns handles to line objects in H.
%
%   Example:
%      Z = eig(randn(20,20));
%      compass(Z)
%
%   See also ROSE, FEATHER, QUIVER.

%   Charles R. Denham, MathWorks 3-20-89
%   Modified, 1-2-92, LS.
%   Modified, 12-12-94, cmt.
%   Copyright 1984-2005 The MathWorks, Inc.

% Parse possible Axes input
narginchk(1,4);
[cax,args,nargs] = axescheck(varargin{:});

xx = [0 1 .8 1 .8].';
yy = [0 0 .08 0 -.08].';
arrow = xx + yy.*sqrt(-1);

if nargs == 2
    x = args{1};
    y = args{2};
    if ischar(y)
        s = y;
        y = imag(x); x = real(x);
    else
        s = [];
    end
elseif nargs == 1
    x = args{1}; 
    s = [];
    y = imag(x); x = real(x);
else % nargs == 3
    [x,y,s] = deal(args{1:3});
end

x = x(:);
y = y(:);
if length(x) ~= length(y)
    error(message('MATLAB:compass:LengthMismatch'));
end

z = (x + y.*sqrt(-1)).';
a = arrow * z;

% Create plot
[th,r] = cart2pol(real(a),imag(a));

if isempty(s),
    h = polar_clean(th,r);
    co = get(cax,'colororder');
    set(h,'color',co(1,:))
else
    h = polar_clean(th,r,s);
end
%if ~isholdon, set(cax,'NextPlot',next); end

% if nargout == 1
%     hc = h;
% end

% Register handles with m-code generator
% if ~isempty(h)
%    mcoderegister('Handles',h,'Target',h(1),'Name','compass');
% end

