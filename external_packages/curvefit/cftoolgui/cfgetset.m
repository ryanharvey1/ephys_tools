function outval = cfgetset(varname,inval)
% CFGETSET is a helper function for CFTOOL

% CFGETSET Get or set a Curve Fitting persistent variable

%   $Revision: 1.11.2.5 $  $Date: 2009/01/23 20:36:42 $
%   Copyright 2001-2008 The MathWorks, Inc.

% Get handle to Curve Fitting figure, usually required later on
oldval = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on');
c = get(0,'Children');
cffig = findobj(c,'flat','Type','figure','Tag','Curve Fitting Figure');
set(0,'ShowHiddenHandles',oldval);

if length(cffig)>1
    cffig = cffig(1);
end

% If asking for the figure itself, handle that now
if isequal(varname,'cffig')
    if isequal(nargin,1)
        outval = cffig;
    elseif ishghandle( inval, 'figure' )
        set(inval,'Tag','Curve Fitting Figure');
    end
    
    % Some things need to persist, so they are root properties
elseif isequal(varname, 'thefitdb')       || isequal(varname, 'thedsdb') || ...
        isequal(varname, 'classinstances') || isequal(varname, 'theoutlierdb')
    
    propname = sprintf( 'CurveFitting_%s', varname );
    
    if nargin==1
        outval = getappdata( 0, propname );
    else
        setappdata( 0, propname, inval );
    end
    
    % For other properties, the figure must not be empty
elseif isempty(cffig)
    if nargout>0
        outval = [];
    end
    
    % If the figure is not empty, and has the appropriate property then set or
    % get it
elseif isprop( cffig, varname )
    if nargin==1
        outval = get( cffig, varname );
    else
        set( cffig, varname, inval );
    end
    
    % Not a builtin property so use appdata
else
    if nargin==1
        outval = getappdata( cffig, lower( varname ) );
    else
        setappdata( cffig, lower( varname ), inval );
    end
end

end

