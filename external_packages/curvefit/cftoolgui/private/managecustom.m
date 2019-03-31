function out = managecustom(action, varargin)
%MANAGECUSTOM A helper function for CFTOOL.

%MANAGECUSTOM Provides access to the custom equations
%
%   OUT = MANAGECUSTOM(ACTION, ...)
%
%   MANAGECUSTOM( 'set', NAME, EQN, OPTS, OLDNAME )
%   MANAGECUSTOM( 'set', NAME, EQN, OPTS, '' )
%
%   The equation, EQN, should be given as a FITTYPE object. The options,
%   OPTS, should be a CURVEFIT.BASEFITOPTIONS object.
%
%   To "edit" an equation, the old name, OLDNAME, for the equation must
%   be supplied. If OLDNAME is empty, then a new entry is made in the
%   library for this equation. For a new equation, the name must not
%   already be in use in the library. If it is in use, then an error
%   will be thrown. Both the "edit" and "new" forms set the dirty flag
%   to true. 
%
%   MANAGECUSTOM( 'delete', NAME ) delete a single equation. Sets the dirty
%   flag to true.
%   MANAGECUSTOM( 'clear') deletes all equations, but doesn't set the dirty
%   flag.
%
%   EQN  = MANAGECUSTOM( 'get', NAME )
%   OPTS = MANAGECUSTOM( 'getopts', NAME )
%   EQNs  = MANAGECUSTOM( 'list' ) is a cell array of all the equations
%   in the library
%
%   STATE = MANAGECUSTOM( 'getdetails', NAME )
%
%   STATE is an 8-by-1 cell array with the following terms,
%     STATE{1} -- formula for the equation, i.e., FORMULA( EQN )
%     STATE{2} -- name of independent variable (char array)
%     STATE{3} -- name of dependent variable (char array)
%     STATE{4} -- coefficient names, i.e., COEFFNAMES( EQN )
%     STATE{5} -- vector of start points
%     STATE{6} -- vector of lower bounds
%     STATE{7} -- vector of upper bounds
%     STATE{8} -- cell-string of linear terms, i.e., linearterms( EQN )
%   If the fit options for this equation do not have start points, lower
%   or upper bounds, then these terms will be empty.

%   $Revision: 1.7.2.8 $  $Date: 2010/10/08 16:36:46 $
%   Copyright 1999-2008 The MathWorks, Inc.

% get (or create) the library
lib = cfgetset( 'customLibrary' );
if isempty( lib )
    lib = i_CreateEmptyLibrary();
end

% Initialize the output
out = '';

% Perform the action
%  -- The "dirty" flag management is done in this switchyard method, NOT
%    in any of the sub-functions below. 
switch action
    case 'set'
        lib = i_Set( lib, varargin{:} );
        cfgetset( 'dirty', true );
    case 'clear'
        lib = i_CreateEmptyLibrary();
    case 'delete'
        lib = i_Delete( lib, varargin{:} );
        cfgetset( 'dirty', true );
    case 'get'
        out = i_GetEquation( lib, varargin{:} );
    case 'getopts'
        out = i_GetOpts( lib, varargin{:} );
    case 'getdetails'
        out = i_GetDetails( lib, varargin{:} );
    case 'list'
        out = lib.names';
end

% save the library
cfgetset( 'customLibrary', lib );

%===============================================================================
function lib = i_Set( lib, name, eq, opts, oldname )
%   MANAGECUSTOM( 'set', NAME, EQN, OPTS, OLDNAME )
%   MANAGECUSTOM( 'set', NAME, EQN, OPTS, '' )

% If the OLDNAME is empty, then this is an 'add'
if isempty( oldname )
    lib = i_AddEquation( lib, name, eq, opts );
else
    % Otherwise it is an edit
    lib = i_EditEquation( lib, oldname, name, eq, opts );
end

%===============================================================================
function lib = i_AddEquation( lib, name, eq, opts )

% Check that an equation with this name doesn't already exist
ind = i_GetEquationIndex( lib, name );
if ~isempty( ind ),
    error(message('curvefit:managecustom:NewEquationNameExists', name));
end

% Set the details of the new equation in the library
lib.names{end+1} = name;
lib.fits{end+1}  = eq;
lib.opts{end+1}  = opts;

%===============================================================================
function lib = i_EditEquation( lib, oldname, newname, eq, opts )

% Check if the old name is in the list
%
% If the old name is not in the list,
%     this is really an "add" so we should do an "add"
%     (this can happen if the equation is deleted while being edited)
%
% Else the old name is on the list
%     If the new name is the same as the old name or is not being used
%         replace the old equation with the new one
%     Otherwise, error

% Check if the old name is in the list
ind = i_GetEquationIndex( lib, oldname );

% If the old name is not in the list,
if isempty( ind )
    % This is really an "add" so we should do an "add"
    lib = i_AddEquation( lib, newname, eq, opts );
    
else % the old name is on the list
    % If the new name is the same as the old name or is not being used
    if isequal( oldname, newname ) || isempty( i_GetEquationIndex( lib, newname ) )
        % Replace the old equation with the new one
        lib.names{ind} = newname;
        lib.fits{ind}  = eq;
        lib.opts{ind}  = opts;
    else
        % The new name is being used ==> error.
        error(message('curvefit:managecustom:EditEquationNameExists', newname));
    end

end

%===============================================================================
function lib = i_CreateEmptyLibrary( )
% MANAGECUSTOM( 'clear')
% Clears the library
lib = struct( ...
    'names', {{}}, ...
    'fits',  {{}}, ...
    'opts',  {{}} );

%===============================================================================
function lib = i_Delete( lib, name )
%   MANAGECUSTOM( 'delete', NAME )

ind = i_GetEquationIndex(lib,name);
lib.names(ind) = [];
lib.fits(ind)  = [];
lib.opts(ind)  = [];

%===============================================================================
function out = i_GetEquation( lib, name )
%   EQN  = MANAGECUSTOM( 'get', NAME )

ind = i_GetEquationIndex(lib,name);
if ~isempty(ind)
    out = lib.fits{ind};
else
    out = [];
end

%===============================================================================
function out = i_GetOpts( lib, name )
%   OPTS = MANAGECUSTOM( 'getopts', NAME )

ind = i_GetEquationIndex(lib,name);
if ~isempty(ind)
    out = copy( lib.opts{ind} );
else
    out = [];
end

%===============================================================================
function out = i_GetDetails( lib, name )
%   STATE = MANAGECUSTOM( 'getdetails', NAME )

ind = i_GetEquationIndex(lib,name);
if ~isempty(ind)
    eq=lib.fits{ind};
    fo=lib.opts{ind};
    % Pull out the first (and only) item in these two, since we
    % can only handle one dependent and one independent variable.
    indname=indepnames(eq);
    depname=dependnames(eq);
    if ~isempty(findprop(fo,'StartPoint'))
        startpoint=fo.StartPoint;
    else
        startpoint=[];
    end
    if ~isempty(findprop(fo,'Lower'))
        lower=fo.Lower;
    else
        lower=[];
    end
    if ~isempty(findprop(fo,'Upper'))
        upper=fo.Upper;
    else
        upper=[];
    end
    out={
        formula(eq)
        indname{1}
        depname{1}
        coeffnames(eq)
        startpoint
        lower
        upper
        linearterms(eq)};
else
    out=[];
end

%=============================================================================== 
function ind = i_GetEquationIndex( lib, name )
ind = strmatch( name, lib.names, 'exact');

% Check that we haven't found more than one equation
if length( ind ) > 1
    % For backward compatibility, just keep the index of the first one
    % found
    ind = ind(1);
    % But still throw a warning because this isn't good.
    warning(message('curvefit:managecustom:RepeatedName', name));
end
