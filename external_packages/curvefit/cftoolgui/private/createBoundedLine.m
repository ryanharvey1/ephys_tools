function h = createBoundedLine(fit, bndsonoff, clev, userargs, varargin)
%CREATEBOUNDEDLINE   Create a "bounded line" object
%
%   H = CREATEBOUNDEDLINE(FIT, BNDSONOFF, CLEV, USERARGS, ...) is a bounded line
%   object.
%
%   FIT is a CFTOOL.FIT object.
%   BNDSONOFF is either the string 'on' or the string 'off'.
%   CLEV is a number between 0 and 100 representing the confidence level to
%   show.
%   USERARGS is a cell array of user arguments. Usually this is {FIT.dataset FIT.dshandle}
%   The variable part of the argument list is parameter-value pairs for
%   properties of a line. These include line specification and parent.

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2010/04/11 20:30:52 $ 

if feature('HGUsingMATLABClasses')
    h = iCreateBoundedLineUsingMATLABClasses(fit,...
        'ShowBounds', bndsonoff, 'ConfLevel', clev, 'UserArgs', userargs,...
        varargin{:} );
else
    h = cftool.boundedline( fit, bndsonoff, clev, '-userargs', userargs,...
        varargin{:} );
end

end

function h = iCreateBoundedLineUsingMATLABClasses( varargin )

if isa(varargin{1},'cftool.fit')
    cftoolFit = varargin{1};
    varargin(1) = [];
else
    cftoolFit = [];
end

h = cftool.BoundedFitLine(varargin{:});
if ~any(strncmpi('parent',varargin,6))
    h.Parent = gca;
end

if ~isempty(cftoolFit)
    h.Fit = cftoolFit;
    h.DFE = cftoolFit.dfe;
    h.String = cftoolFit.name;
end


end