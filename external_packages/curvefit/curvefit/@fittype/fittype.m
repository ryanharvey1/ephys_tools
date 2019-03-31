classdef fittype
    %FITTYPE   Fittype for curve and surface fitting
    %
    %   A FITTYPE encapsulates information describing a model. To create a
    %   fit, you need data, a FITTYPE, and (optionally) fit options and an
    %   exclusion rule.
    %
    %
    %   LIBRARY MODELS
    %
    %   FITTYPE(LIBNAME) constructs a FITTYPE for the library model
    %   specified by LIBNAME.
    %
    %   Choices for LIBNAME include:
    %
    %       LIBNAME           DESCRIPTION
    %       'poly1'           Linear polynomial curve
    %       'poly11'          Linear polynomial surface
    %       'poly2'           Quadratic polynomial curve
    %       'linearinterp'    Piecewise linear interpolation
    %       'cubicinterp'     Piecewise cubic interpolation
    %       'smoothingspline' Smoothing spline (curve)
    %       'lowess'          Local linear regression (surface)
    %
    %   or any of the names of library models described in CFLIBHELP (type
    %   CFLIBHELP to see the names and descriptions of library models).
    %
    %
    %   CUSTOM MODELS
    %
    %   FITTYPE(EXPR) constructs a FITTYPE from the MATLAB expression
    %   contained in the string, cell array or anonymous function EXPR.
    %
    %   The FITTYPE automatically determines input arguments by searching
    %   EXPR for variable names (see SYMVAR). In this case, the FITTYPE
    %   assumes 'x' is the independent variable, 'y' is the dependent
    %   variable, and all other variables are coefficients of the model. If
    %   no variable exists, 'x' is used.
    %
    %   All coefficients must be scalars. You cannot use the following
    %   coefficient names in the expression string EXPR: i, j, pi, inf,
    %   nan, eps.
    %
    %   If EXPR is a string or anonymous function, then the toolbox uses a
    %   nonlinear fitting algorithm to fit the model to data. To use a
    %   linear fitting algorithm, use a cell array of terms.
    %
    %
    %   ANONYMOUS FUNCTIONS
    %
    %   If EXPR is an anonymous function, then the order of inputs must be
    %   correct. The input order enables the FITTYPE class to determine
    %   which inputs are coefficients to estimate, problem-dependent
    %   parameters and independent variables. The order of the input
    %   arguments to the anonymous function must be:
    %
    %   EXPR = @(<coefficients>, <problem parameters>, <x>, <y>) expression
    %
    %   There must be at least one coefficient. The problem parameters and
    %   y are optional. The last arguments, x and y, represent the
    %   independent variables: just x for curves, but x and y for surfaces.
    %   If you don't want to use x and/or y as the names of the independent
    %   variables, then you can specify different names by using the
    %   'independent' property name/value pair. However, whatever name or
    %   names you choose, these arguments must be the last arguments to the
    %   anonymous function.
    %
    %   Anonymous functions make it easier to pass other data into the
    %   fittype and fit functions. For example, to create a fittype using
    %   an anonymous function and variables xs and ys from the workspace:
    %
    %       % Variables in workspace
    %       xs = (0:0.1:1).';
    %       ys = [0; 0; 0.04; 0.1; 0.2; 0.5; 0.8; 0.9; 0.96; 1; 1];
    %       % Create fittype
    %       ft = fittype( @(b, h, x) interp1( xs, b+h*ys, x, 'pchip' ) )
    %       % Load some data
    %       xdata = [0.012;0.054;0.13;0.16;0.31;0.34;0.47;0.53;0.53;...
    %           0.57;0.78;0.79;0.93];
    %       ydata = [0.78;0.87;1;1.1;0.96;0.88;0.56;0.5;0.5;0.5;0.63;...
    %           0.62;0.39];
    %       % Fit the curve to the data
    %       f = fit( xdata, ydata, ft, 'Start', [0, 1] )
    %
    %
    %   LINEAR MODELS
    %
    %   To use a linear fitting algorithm specify EXPR as a cell array of
    %   terms. That is, to specify a linear model of the following form:
    %
    %       coeff1 * term1 + coeff2 * term2 + coeff3 * term3 + ...
    %
    %   (where no coefficient appears within any of term1, term2, etc) use
    %   a cell array where each term, without coefficients, is specified in
    %   a cell of EXPR, as follows:
    %
    %       EXPR = {'term1', 'term2', 'term3', ... }
    %
    %   For example, the model
    %
    %       a*x + b*sin(x) + c
    %
    %   is linear in 'a', 'b' and 'c'. It has three terms 'x', 'sin(x)' and
    %   '1' (since c=c*1) and so EXPR is
    %
    %       EXPR = {'x','sin(x)','1'}
    %
    %
    %   ADDITIONAL PROPERTIES
    %
    %   FITTYPE(EXPR,PROP1,VALUE1,PROP2,VALUE2,....) uses the property
    %   name/value pairs PROP1-VALUE1, PROP2-VALUE2 to specify property
    %   values other than the default values.
    %
    %   PROPERTY         DESCRIPTION
    %   'independent'    Specifies the independent variable name
    %   'dependent'      Specifies the dependent variable name
    %   'coefficients'   Specifies the coefficient names (in a cell array
    %                    if there are two or more). Note excluded names 
    %                    above.
    %   'problem'        Specifies the problem-dependent (constant) names
    %                    (in a cell array if there are two or more)
    %   'options'        Specifies the default 'FITOPTIONS' for this 
    %                    equation
    %
    %   Defaults: The independent variable is x.
    %             The dependent variable is y.
    %             There are no problem dependent variables.
    %             Everything else is a coefficient name.
    %
    %   Multi-character symbol names may be used.
    %
    %
    %   EXAMPLES
    %
    %      g = fittype('a*x^2+b*x+c')
    %      g = fittype('a*x^2+b*x+c','coeff',{'a','b','c'})
    %      g = fittype('a*time^2+b*time+c','indep','time')
    %      g = fittype('a*time^2+b*time+c','indep','time','depen','height')
    %      g = fittype('a*x+n*b','problem','n')
    %      g = fittype({'cos(x)','1'})                            % linear
    %      g = fittype({'cos(x)','1'}, 'coefficients', {'a','b'}) % linear
    %      g = fittype( @(a,b,c,x) a*x.^2+b*x+c )
    %      g = fittype( @(a,b,c,d,x,y) a*x.^2+b*x+c*exp(-(y-d).^2), ...
    %           'independent', {'x', 'y'}, ...
    %           'dependent', 'z' ); % for fitting surfaces
    %
    %   See also cflibhelp, fit, fitoptions.
    
    %   Copyright 1999-2011 The MathWorks, Inc.
    %   $Revision: 1.25.2.25 $  $Date: 2011/05/09 00:39:24 $
    
    properties (Access=private)
        % fType   string of type information
        fType = '';
        
        % fTypename   string description for display of fittype
        fTypename = 'Model';
        
        % fCategory   string of category information
        fCategory = '';
        
        % defn   string containing function definition for display
        defn = '';
        
        % fFeval   flag indicating if model is evaluated using feval or eval
        fFeval = 0;
        
        % expr   string expression to eval or function (handle) to feval
        expr = '';
        
        % Adefn    cell array of strings to display A matrix for linear custom
        % equations
        Adefn = {};
        
        % Aexpr   cell array of strings to form A matrix for linear custom
        % equations
        Aexpr = {};
        
        % linear   flag indicating if model is linear
        linear = 0;
        
        % derexpr   string expression or function (handle) of derivative
        derexpr = '';
        
        % intexpr   string expression or function (handle) of integral
        intexpr = '';
        
        % args   string matrix containing formal parameter names
        args = '';
        
        % isEmpty   flag indicating fittype called with no arguments
        isEmpty = 1;
        
        % numArgs   number of formal parameters
        numArgs = 0;
        
        % numCoeffs   number of coefficients (unknowns)
        numCoeffs = 0;
        
        % assignData   eval string to assign inputs to independent (data)
        % variable
        assignCoeff= '';
        
        % assignCoeff   eval string to assign inputs to coefficients
        assignData= '';
        
        % assignProb   eval string to assign inputs to problem dependent
        % parameters
        assignProb= '';
        
        % indep   string matrix of independent parameter names (one name per
        % row)
        indep = '';
        
        % depen   string matrix of dependent parameter name
        depen = '';
        
        % coeff   string matrix of coefficient parameter names (one name per row)
        coeff = '';
        
        % prob   string matrix of problem parameter names (one name per row)
        prob = '';
        
        % fConstants   cell array of "hidden" constants that are passed to
        % feval'd library function after problem parameters, but are not counted
        % in .numArgs (hidden from user).
        fConstants = {};
        
        % nonlinearcoeff   index of nonlinear coefficients in separable
        % equations
        fNonlinearcoeffs = [];
        
        % fFitoptions   default options to use in fitting
        fFitoptions = '';
        
        % fStartpt   function (handle) of starting point computer
        fStartpt = '';
        
        % version   object version
        version = 2.0;
    end
    
    methods(Access = 'protected', Hidden)
        varargout = evaluate(varargin)
        line = fcnstring( obj, variable, numargs, arglist )
        line = nonParametricFcnString( obj, variable, arglist, coefficient )
        obj_out = copyFittypeProperties( obj_out, obj_in )
    end
    
    methods(Static, Hidden)
        obj = loadobj(obj)
    end
    
    methods(Hidden)
        c = cat(varargin)
        str = char(obj)
        obj = clearhandles(obj)
        constants = constants(fun)
        derexpr = derivexpr(fun)
        disp(obj)
        c = exist(obj)
        args = fevalexpr(fun)
        FITTYPE_OBJ_A_ = getcoeffmatrix(FITTYPE_OBJ_,varargin)
        c = horzcat(varargin)
        intexpr = integexpr(fun)
        c = isempty(obj)
        cella = linearexprs(model)
        cella = linearterms(model)
        y = nargin(obj)
        y = nargout(obj)
        cella = nonlinearcoeffs(model)
        n = numindep(fun)
        type = prettyname(fittype)
        obj = saveobj(obj)
        num = startpt(model)
        c = subsasgn(FITTYPE_OBJ_, varargin)
        FITTYPE_OUT_ = subsref(FITTYPE_OBJ_, FITTYPE_SUBS_)
        vars = symvar(s)
        c = vertcat(varargin)
    end
    
    methods
        function obj = fittype(varargin)
            %FITTYPE  Construct a fit type object.
            %
            %   See also: fittype.
            
            %  Notes:
            %  1. Coefficients must all be scalars:
            %     g = fittype('w(1)*x^2+w(2)', 'coeff','w') will error because in SYMVAR, 'w'
            %     won't be detected as a variable (it looks like a function).
            %  2. Linear custom equations are created by passing in a cell array of
            %     expressions as the first argument.  In this case, the 'coefficients'
            %     param can be provided with a list of coefficient names in the order
            %     corresponding to the expressions in the cell array.  If there is a
            %     constant term, use '1' as the corresponding expression in the cell array.
            %     Example:
            %       g = fittype({'cos(x)','1'})
            %       or
            %       g = fittype({'cos(x)','1'}, 'coefficients', {'a','b'})
            %     sets up the custom equation 'a*cos(x) + b', and will treat it as a
            %     linear equation.
            %  3. Library functions can have "hidden" constants that can be set.
            %     See .fConstants below.
            %     Example:
            %       g = fittype('rat34')
            %     will set up function handles to library functions that need to get the 3,
            %     4 as constants to compute.
            
            % Some potential early exit
            if nargin==0
                % return default object
                return
            elseif isequal(class(varargin{1}),'fittype') && nargin==1
                % if fittype object, just return
                obj = varargin{1};
                return
            elseif isa(varargin{1},'fittype') && nargin==1
                % else if subclass of fittype object, return fittype
                obj = copyFittypeProperties( obj, varargin{1} );
                return
            end
            
            % First input argument determines the type of fittype to create
            if isidentifier(varargin{1})
                % look for spline, interpolant or library function
                obj = iCreateFromLibrary( obj, varargin{:} );
                
            elseif isempty(varargin{1})
                % fittype('',...)
                error(message('curvefit:fittype:emptyExpression'));
                
            elseif iscellstr(varargin{1}) || ischar(varargin{1})
                % turn expression into a custom model: nonlinear or linear
                obj = iCreateCustomFittype( obj, varargin{:} );
                
            elseif iIsAnonymousFunction( varargin{1} )
                % Turn anonymous function into a custom fittype
                obj = iCreateFromAnonymousFunction( obj, varargin{:} );
                
            else
                % otherwise first argument is something that we can't deal
                % with.
                error(message('curvefit:fittype:firstArgMustBeStringOrCell'));
            end
            
            % Re-order the arguments to the fittype so that we know how to
            % evaluate
            obj = iSetArguments( obj );
            
            % For custom models, we need to check that the parameters are all
            % accounted for and that we can actually evaluate the model.
            if isequal( 'custom', category( obj ) )
                iTestCustomModelParameters( obj );
                iTestCustomModelEvaluation( obj );
            end
        end
    end
    
end

%------------------------------------------
% Private methods that should only be used during construction
%------------------------------------------
function obj = iCreateFromLibrary( obj, varargin )
% iCreateFromLibrary -- create fittype from library model
p = inputParser;
p.CaseSensitive = false;
p.KeepUnmatched = true;
p.FunctionName = 'fittype';
p.addParamValue( 'numindep', [], @(n) isnumeric( n ) && (isscalar( n ) || isempty( n )) );
p.parse( varargin{2:end} );

% Check for unmatched parameter-value pairs
if ~isempty(fieldnames( p.Unmatched ) )
    error(message('curvefit:fittype:TooManyInputsForLibraryModel'));
end

nindep = p.Results.numindep;

[fittypes,fitcateg,fitnames,allowedNIndep] = getfittypes;
% Start by looking for an exact match and if that fails, look
% for a unique partial match
ind = find( strcmpi( varargin{1}, cellstr( fittypes ) ) );
if isempty( ind )
    ind = find( strncmpi( varargin{1}, cellstr( fittypes ), length( varargin{1} ) ) );
end
if length(ind) > 1
    error(message('curvefit:fittype:typeNotUnique', varargin{ 1 }));
end
if ~isempty(ind) && (2 <= ind)
    obj.isEmpty = 0;
    obj.fCategory = deblank(fitcateg(ind,:));
    obj.fType = deblank(fittypes(ind,:));
    obj.fTypename = deblank(fitnames(ind,:));
    obj = iIndepAndDepenFromNumIndep( obj, nindep, allowedNIndep{ind} );
    
    % get defn, expr, derexpr, intexpr, feval, coeff, or error out
    obj = liblookup( obj.fType, obj );
else
    error(message('curvefit:fittype:fcnNotFound', varargin{ 1 }));
end
end

%------------------------------------------
function obj = iCreateCustomFittype( obj, varargin )
% iCreateCustomFittype -- create fittype for custom model.
obj.fCategory = 'custom';
% Get input expression and formal parameters
obj.isEmpty = 0;
% we assume x and y unless other parameters change this
obj.indep = 'x';
obj.depen = 'y';

% Parse param-value pairs
obj = iParseParameters(obj,varargin(2:end));

if ischar(varargin{1}) % nonlinear custom equation
    obj.fType = 'customnonlinear';
    obj.defn = strtrim(varargin{1});

    if isempty(obj.fFitoptions)
        obj.fFitoptions = fitoptions('method','nonlinearleastsquares');
    end
    
else % cell array: linear custom equation
    obj.fType = 'customlinear';
    obj.linear = 1;
    
    % Derive the definition and expression for A from the first input argument
    [obj.Adefn, obj.Aexpr] = iFilterAdefn( varargin{1} );
    
    % Generate coefficients is none supplied
    if isempty(obj.coeff)
        obj.coeff = iDefaultLinearCoefficients( obj );
        
        % ... otherwise check that the user has given the right number of
        % coefficients
    elseif size(obj.coeff,1) ~= length(obj.Aexpr)
        error(message('curvefit:fittype:numCoeffsAndTermsMustMatch'));
    end
    
    % Derive definition of model from terms and coefficients
    obj.defn = iCustomLinearDefinition( obj.Adefn, obj.coeff );
    
    if isempty(obj.fFitoptions)
        obj.fFitoptions = fitoptions('method','linearleastsquares');
    end
end

obj.expr = vectorize(obj.defn);
obj.args = char( symvar( obj.expr ) );
% Take care of the case:  FITTYPE('2').
if isempty(obj.args)
    obj.args = 'x';
end

% Use heuristics to get formal parameters
if isempty(obj.coeff)
    obj = iDeduceCoefficients(obj);
end
obj.numCoeffs = size(obj.coeff,1);

end

%------------------------------------------
function obj = iCreateFromAnonymousFunction( obj, varargin )
% iCreateFromAnonymousFunction -- Create a fittype from an Anonymous
% Function.

% The first argument is the function handle to base the fittype on
theFcn = varargin{1};

% The anonymous function is not allowed to use varargin
if nargin( theFcn ) < 0
    error(message('curvefit:fittype:AnonFcnVarargin'));
end

% This is not an empty fittpye
obj.isEmpty = 0;

% We are going to treat an anonymous function as a "custom non-linear"
% model.
obj.fType = 'customnonlinear';
obj.fCategory = 'custom';

% By default, the independent variable is 'x'
obj.indep = 'x';
% By default, the dependent variable is 'y'
obj.depen = 'y';

% Parse param-value pairs
obj = iParseParameters( obj, varargin(2:end) );

% If no coefficients are given,
if isempty( obj.coeff )
    % ... then we need to deduce them from the argument list
    obj.args = char( iArgumentsFromAnonymousFunction( theFcn ) );
    obj = iDeduceCoefficients(obj);
end
obj.numCoeffs = size( obj.coeff, 1 );

% Set the formula (or definition) 
obj.defn = iFormulaFromAnonymousFunction( theFcn );

% To set up evaluation, 
% ... need to set the fFeval flag to true, i.e., use "feval"
obj.fFeval = true;
% ... and set the expression, expr, to be the function handle
obj.expr = theFcn;

% Specify the fit options to use
if isempty( obj.fFitoptions )
    obj.fFitoptions = fitoptions( 'method', 'nonlinearleastsquares' );
end

% Check order of arguments
iTestAnonymousFunctionArgumentOrder( obj, theFcn )
end

%------------------------------------------
function obj = iParseParameters(obj,arglist)
%iParseParameters --  parse parameter value pairs.

% Check the we have parameter-value PAIRS.
n = length(arglist);
if mod(n,2)
    error(message('curvefit:fittype:parameterNamesAndValuesNotInPairs'));
end

% List of all valid parameter names
ALL_PARAMETERS = {
    'independent'
    'dependent'
    'coefficients'
    'numindep'
    'problem'
    'options'
    };

% Process each parameter in turn
for i = 1:2:n
    param = arglist{i};
    value = arglist{i+1};

    % Look for unique partial match on the parameter name
    tf = strncmpi( param, ALL_PARAMETERS, length( param ) );
    if nnz( tf ) ~= 1
        iInvalidParameterError( param );
    end
    % Use the full name of the parameter
    param = ALL_PARAMETERS{tf};
    
    switch param
        case 'independent'
            obj.indep = char(value); 
        case 'dependent'
            obj.depen = char(value); 
        case 'coefficients'
            obj.coeff = char(value); 
        case 'numindep'
            % For the moment, we assume that this is a valid number of
            % independent variables.
            obj = iIndepAndDepenFromNumIndep( obj, value );
        case 'problem'
            obj.prob = char(value); 
        case 'options'
            obj.fFitoptions = copy( value );
        otherwise
            % Shouldn't hit this case.
            iInvalidParameterError( param );
    end % switch
end

end

%------------------------------------------
function obj = iIndepAndDepenFromNumIndep( obj, nindep, allowedNIndep )
% iIndepAndDepenFromNumIndep -- Set names of independent and dependent variables
% based on the number of independent variables
%
% Also check that the number of independent variables requested is allowed for
% this model type.
if nargin < 3
    allowedNIndep = [1, 2];
end

if isempty( nindep )
    % By default take the first allow number
    nindep = allowedNIndep(1);
    
elseif ~ismember( nindep, allowedNIndep )
    iInvalidNumIndepError( allowedNIndep, obj.fTypename );
end

switch nindep
    case 1
        obj.indep = 'x';
        obj.depen = 'y';
    case 2
        obj.indep = ['x';'y'];
        obj.depen = 'z';
    otherwise
        ME = MException( 'curvefit:fittype:InvalidNindep', ...
            'Number of independent variables, NINDEP, must be either 1 or 2.' );
        throwAsCaller( ME );
end
end

%------------------------------------------
function obj = iDeduceCoefficients(obj)
%iDeduceCoefficients Figure out coeff names from obj.args and obj.indep

% Start by assuming that all arguments are coefficients
index = true(1,size(obj.args,1));

% Look for elements of the argument array that are independent variables
for j = 1:size(obj.indep,1)
    i = strcmp( strtrim( obj.indep(j,:) ), cellstr( obj.args ) );
    % If this name isn't found, then error
    if ~any( i )
        error(message('curvefit:fittype:missingIndVar', obj.indep( j, : ), obj.indep( j, : )));
    end
    % ... else remove that name from the list of coefficients
    index(i) = false;

end

% Looks for elements of the argument array that are problem parameters
for j = 1:size(obj.prob,1)
    i = strcmp( strtrim( obj.prob(j,:) ), cellstr( obj.args ) );
    % If this name isn't found, then error
    if ~any( i )
        error(message('curvefit:fittype:invalidProbParam'))
    end
    % ... else remove that name from the list of coefficients
    index(i) = false;
end

% If there are no (true) entries in the list of coefficients, then there are no
% coefficients in the expression and we have to throw an error.
if ~any( index )
    error(message('curvefit:fittype:noCoeffs'));
end

% Set the coefficient property with the names of the coefficients that we found.
obj.coeff = obj.args(index,:);

end

%------------------------------------------
function coeff = iDefaultLinearCoefficients( obj )

numcoeff = length(obj.Aexpr);

% The default list of coefficients is the whole alphabet minus i and j
coeffnames = {'a';'b';'c';'d';'e';'f';'g';'h';'k';'l';'m'; ...
    'n';'o';'p';'q';'r';'s';'t';'u';'v';'w';'x';'y';'z'};

% We now look for any names already used that clash with these coefficients
usedNames = cellstr( [obj.indep; obj.depen; obj.prob] );
clashes = ismember( coeffnames, usedNames );
% ... and remove them from the list
coeffnames(clashes) = [];

% If there are insufficient default coefficient names for the number of
% terms provided by the user...
if numcoeff > length( coeffnames )
    % ... then throw an error
    error(message('curvefit:fittype:tooManyCoeffNames'))
end
% Return one coefficient per term in the model.
% (Return variable needs to be a char-array.)
coeff = char( coeffnames(1:numcoeff,:) );

end

%------------------------------------------
function obj = iSetArguments( obj )
% iSetArguments -- Set the "arguments" of the FITTYPE.
%
% This function sets the "args" and "numArgs" properties of the FITTYPE based on
% the coefficients, problem parameters and independent variables.
%
% It also sets up the mappings (assignCoeff, assignProb and assignData) between
% the argument names and inputs to the feval method.

% NOTE: For custom models, we have already set the args property, but for
% library models it needs to be set.

if isempty( obj.prob )
    obj.args = char(obj.coeff, obj.indep); 
else
    obj.args = char(obj.coeff, obj.prob, obj.indep); 
end
obj.numArgs = size(obj.args,1);

% Mappings for coefficients
numCoeff = iNumArguments( obj.coeff );
for j = 1:numCoeff
    obj.assignCoeff = sprintf('%s %s = FITTYPE_INPUTS_{%d};', ...
        obj.assignCoeff, deblank(obj.coeff(j,:)), j);
end

% Mappings for problem parameters
numProb = iNumArguments( obj.prob );
for m = 1:numProb
    obj.assignProb = sprintf('%s %s = FITTYPE_INPUTS_{%d};', ...
        obj.assignProb, deblank(obj.prob(m,:)), numCoeff+m);
end

% Mappings for independent variables
numIndep = iNumArguments( obj.indep );
for k = 1:numIndep
    obj.assignData = sprintf('%s %s = FITTYPE_INPUTS_{%d};',...
        obj.assignData, deblank(obj.indep(k,:)), numCoeff+numProb+k);
end

end

%------------------------------------------
function iTestCustomModelEvaluation( obj )

try
    xtemp = arbitrary( 2 ).'; % create data
    probparams = num2cell( arbitrary( numargs(obj)-obj.numCoeffs-1 ) );
    
    if islinear(obj) % Also eval to get coefficient matrix
        ignore = getcoeffmatrix(obj,probparams{:},xtemp); %#ok<NASGU>
    end
    
    if ~isempty(findprop(obj.fFitoptions, 'StartPoint')) && ...
            isequal(length(obj.fFitoptions.StartPoint),obj.numCoeffs)
        vals = num2cell(obj.fFitoptions.StartPoint);
    else
        vals = num2cell( arbitrary( obj.numCoeffs ) );
    end
    
    % try to evaluate the model
    ignore = feval(obj,vals{:},probparams{:},xtemp); %#ok<NASGU>
    
catch caughtError
    newError = MException( 'curvefit:fittype:invalidExpression', ...
        'Expression %s is not a valid MATLAB expression,\n has non-scalar coefficients, or cannot be evaluated:\n%s',...
        obj.defn, caughtError.message );
    addCause( newError, caughtError );
    throwAsCaller( newError );
end
end

%------------------------------------------
function iTestCustomModelParameters( obj )
% Check that all coeff are in an expression and that extra
% variables are not coeffs

% Get a list of all the variables in the model
vars = symvar( obj.defn );

% obj.coeff, obj.prob and obj.indep are char arrays
% possibly with trailing blanks. Convert to cell arrays
% without blanks to make exact comparison possible.
if ~isempty(obj.indep)
    indepcell = deblank(cellstr(obj.indep));
else
    indepcell = {};
end
if ~isempty(obj.coeff)
    coeffcell = deblank(cellstr(obj.coeff));
else
    coeffcell = {};
end
if ~isempty(obj.prob)
    probcell = deblank(cellstr(obj.prob));
else
    probcell = {};
end

for i=1:length(indepcell)
    if ~ismember( indepcell{i}, vars )
        error(message('curvefit:fittype:noIndependentVar', indepcell{ i }));
    end
end
for i=1:length(coeffcell)
    if ~ismember( coeffcell{i}, vars )
        error(message('curvefit:fittype:noCoeff', coeffcell{ i }));
    end
end
for i=1:length(probcell)
    if ~ismember( probcell{i}, vars )
        error(message('curvefit:fittype:noProbParam', probcell{ i }));
    end
    % compare prob and coeff
    if ismember(probcell{i},coeffcell)
        error(message('curvefit:fittype:SameNameCoeffAndProb', probcell{ i }));
    end
end
for i=1:length( indepcell )
    if ismember( indepcell{i}, coeffcell )
        error(message('curvefit:fittype:SameNameCoeffAndIndVar', indepcell{ 1 }));
    end
    if ismember( indepcell{i}, probcell )
        error(message('curvefit:fittype:SameNameProbParamAndIndVar', indepcell{ i }));
    end
end
if ismember(obj.depen,coeffcell)
    error(message('curvefit:fittype:SameNameCoeffAndDepenVar', obj.depen));
end
if ismember(obj.depen,probcell)
    error(message('curvefit:fittype:SameNameProbParamAndDepenVar', obj.depen));
end
if ismember(obj.depen,cellstr( obj.indep ))
    error(message('curvefit:fittype:SameNameIndAndDepenVars', obj.indep));
end

% Check for duplicate argument names
if numel( coeffcell ) ~= numel( unique( coeffcell ) )
    error(message('curvefit:fittype:duplicateCoefficients'));
end
if numel( probcell ) ~= numel( unique( probcell ) )
    error(message('curvefit:fittype:duplicateProblemParameters'));
end
if numel( indepcell ) ~= numel( unique( indepcell ) )
    error(message('curvefit:fittype:duplicateIndependent'));
end
% Not only are you not allowed duplicate dependent variable names, you are only
% allowed one dependent variable
if size( obj.depen, 1 ) ~= 1
    error(message('curvefit:fittype:multipleDependent'));
end

% Checking if coeffs appear in linear terms for linear custom
if obj.linear
    linearvars = {};
    for i=1:length(obj.Adefn)
        linearvars = vertcat( linearvars, symvar( obj.Adefn{i} ) ); %#ok<AGROW>
    end
    for i=1:size(obj.coeff,1)
        if ismember( obj.coeff(i,:), linearvars )
            error(message('curvefit:fittype:linearTermContainsCoeff', obj.coeff( i, : )));
        end
    end
end
end

%------------------------------------------
function iTestAnonymousFunctionArgumentOrder( obj, theFcn )
% iTestAnonymousFunctionArgumentOrder -- check the order of the inputs to
% an anonymous function:
%
% 1. The coefficients have to come before the independent variables
%
% 2. The coefficients have to come before any problem parameters
%
% 3. Any problem parameters have to come before the independent variables
%
% 4. The independent variables have to appear in the same order in the
% anonymous functions as the do in the parameter-value pair listing
%
% 5. Any problem parameters must also appear in the same order as in the
% parameter-value pair listing
%
% 6. The coefficients have to appear in the same order in the anonymous
% functions as the do in the parameter-value pair listing
%
% Order:
%   <coefficients> <problem parameters> <independent variables>

arguments = iArgumentsFromAnonymousFunction( theFcn );
[~, coefficientLocations] = ismember( obj.coeff, arguments );
[~, independentLocations] = ismember( obj.indep, arguments );

% 1. The coefficients have to come before the independent variables
if max( coefficientLocations ) > min( independentLocations )
    error(message('curvefit:fittype:AnonFcn:IndepBeforeCoeffs'));
end

% 6. The coefficients have to appear in the same order in the anonymous
% functions as the do in the parameter-value pair listing
if ~issorted( coefficientLocations )
    error(message('curvefit:fittype:AnonFcn:WrongCoefficientsOrder'));
end

% 4. The independent variables have to appear in the same order in the anonymous
% functions as the do in the parameter-value pair listing
if ~issorted( independentLocations )
    error(message('curvefit:fittype:AnonFcn:WrongIndependentVariablesOrder'));
end

if ~isempty( obj.prob )
    [~, problemLocations] = ismember( obj.prob,  arguments );

    % 2. Any  coefficients have to come before the problem parameters
    if max( coefficientLocations ) > min( problemLocations )
        error(message('curvefit:fittype:AnonFcn:ProblemBeforeCoeffs'));
    end
    
    % 3. Any problem parameters have to come before the independent variables
    if max( problemLocations ) > min( independentLocations )
        error(message('curvefit:fittype:AnonFcn:IndepBeforeProblem'));
    end
    
    % 5. Any problem parameters must also appear in the same order as in the
    % parameter-value pair listing 
    if ~issorted( problemLocations )
        error(message('curvefit:fittype:AnonFcn:WrongProblemParametersOrder'));
    end
end
end

%------------------------------------------
% Helper functions
%------------------------------------------
function s1 = strtrim(s)
%STRTRIM Trim spaces from string.

if isempty(s)
    s1 = s;
else
    % remove leading and trailing blanks (including nulls)
    c = find(s ~= ' ' & s ~= 0);
    s1 = s(min(c):max(c));
end
end

%-------------------------------------------
function tf = isidentifier(str)

tf = 0;
if ~ischar(str)
    return;
end

if ~isempty(str)
    first = str(1);
    if (isletter(first))
        letters = isletter(str);
        numerals = (48 <= str) & (str <= 57);
        underscore = (95 == str);
        if (all(letters | numerals | underscore))
            tf = 1;
        end
    end
end

tf = logical(tf);

end

%------------------------------------------
function iInvalidNumIndepError( allowedNIndep, fitname )
% iInvalidNumIndepError -- Throw an appropriate error for an invalid number
% of independent variables.

% The error messages will contain a reference to CFLIBHELP
if feature( 'hotlinks' ) && ~isdeployed
    linkToCflibhelp = sprintf( '<a href="matlab:cflibhelp">cflibhelp</a>' );
else
    linkToCflibhelp = 'cflibhelp';
end

% The expected values for "allowedNIndep" are 1, 2 and [1, 2]. There is a
% different exception for each of these cases.
if isequal( allowedNIndep, 1 )
    theException = MException( 'curvefit:fittype:InvalidNumIndep1', ...
        sprintf( 'The library function ''%s'' is for fitting curves and requires one independent variable, that is, X must have one column.\nSee %s for information on all library models.', ...
        fitname, linkToCflibhelp ) );

elseif isequal( allowedNIndep, 2 )
    theException = MException( 'curvefit:fittype:InvalidNumIndep2', ...
        sprintf( 'The library function ''%s'' is for fitting surfaces and requires two independent variables, that is, X must have two columns.\nSee %s for information on all library models.', ...
        fitname, linkToCflibhelp ) );
    
else 
    % Assume that allowedNIndep holds two valid values
    theException = MException( 'curvefit:fittype:InvalidNumIndep12', ...
        sprintf( 'The library function ''%s'' is for fitting models with %d or %d independent variables, that is, X must have %d or %d columns.\nSee %s for information on all library models.', ...
        fitname, allowedNIndep(1), allowedNIndep(2), allowedNIndep(1), allowedNIndep(2), linkToCflibhelp ) );
end
throwAsCaller( theException );
end

%------------------------------------------
function r = arbitrary( k )
% ARBITRARY( K ) is a row vector of K arbitrary values. This can be
% thought of as a substitute for RAND( 1, K ) when the "random" aspect
% is not required.
r = (1:k)/(k+1);
end

%------------------------------------------
function num = iNumArguments( arguments )
% iNumArguments -- Count the number of "arguments" is char array.
%
% There is one argument per row, unless the char array is empty in which case it
% will have one row but no arguments.
if isempty( arguments )
    num = 0;
else
    num = size( arguments, 1 );
end
end

%------------------------------------------
function [Adefn, Aexpr] = iFilterAdefn( Adefn )

Aexpr = cell( size( Adefn ) );

for i=1:length( Adefn )
    if ~isempty( Adefn{i} )
        if ~isempty( [strfind( Adefn{i},'+'),strfind( Adefn{i},'-')] )
            % + or - in expression: add parenthesis
            Adefn{i} = ['(',Adefn{i},')'];
        end
        Aexpr{i} = vectorize(Adefn{i});
    else
        error(message('curvefit:fittype:emptyLinearTerm'));
    end
end

end

%------------------------------------------
function definition = iCustomLinearDefinition( Adefn, coeff )

definition = '';
for i=1:length(Adefn)
    if isequal(Adefn{i},'1')
        definition = [definition, deblank(coeff(i,:)),' + ']; %#ok<AGROW>
    else
        definition = [definition, deblank(coeff(i,:)),'*',Adefn{i}, ' + ']; %#ok<AGROW>
    end
end
definition = definition(1:end-3); % trim off that last +

end

%------------------------------------------
function arguments = iArgumentsFromAnonymousFunction( theFcn )
% iArgumentsFromAnonymousFunction -- A cell-string of the input arguments
% to an anonymous function.
%
% The order of the strings in the cell array is the same as the order of
% the arguments to the anonymous function.
%
% If you call |func2str| on an anonymous function then the resulting string
% like this:
%
%     <some text>@(<arg1>,<arg2>,...,<argN>)<text of function>
%f-
% We want to get the names 'arg1', 'arg2', ..., 'argN' and return them as a
% cell-string.
%
% Also, we can know the number of arguments, N, ahead of time by calling
% "nargin" of the function handle.
numArguments = nargin( theFcn );

% Get the string description of the function
functionString = func2str( theFcn );

% Allocate space for the cell-string
arguments = cell( numArguments, 1 );

% The plan is to move a pair of indices along the "function string" field
% looking for the commas. The names we want will be between these indices.
% 
% We know from the form of the string that the first name starts two
% characters after the "@"
indexOfAtSign = find( functionString == '@', 1, 'first' );
ai = indexOfAtSign + 2;
% Therefore the first comma must be no sooner that the fourth character
bi = 4;
% When we start we have found no arguments
numFound = 0;
% We will keep looping until we have found all the arguments we expect
while numFound < numArguments
    % If we have found the end of a argument name
    if functionString(bi) == ',' || functionString(bi) == ')'
        % then increment the "numFound" counter 
        numFound = numFound+1;
        % ... store the name
        arguments{numFound} = functionString(ai:(bi-1));
        % and increment the start index
        ai = bi+1;
        % Since the end must be beyond the start, we set the end index
        % beyond the start index.
        bi = ai+1;
    else
        % Otherwise increment the end index
        bi = bi+1;
    end
end
end

%------------------------------------------
function formula = iFormulaFromAnonymousFunction( theFcn )
% iFormulaFromAnonymousFunction -- A string that describes an anonymous
% function.
%
% We look at the string returned by the |func2str| function applied to the
% anonymous function. Everything after the "@(...)" is the formula

% Get the string description of the function
functionString = func2str( theFcn );

% Find the index of the first closing bracket
idx = find( functionString == ')', 1, 'first' );

% The formula is everything from this point on
formula = functionString((idx+1):end);
end

%------------------------------------------
function iInvalidParameterError( parameter )
if ischar( parameter )
    error(message('curvefit:fittype:invalidParam', parameter));
else
    error(message('curvefit:fittype:parameterNameNotAString'));
end
end

%------------------------------------------
function tf = iIsAnonymousFunction( obj )
% iIsAnonymousFunction -- True if object is an anonymous function
%
% An object is an anonymous function if it is a function handle and the
% first two characters of the string form are '@('
tf = isa( obj, 'function_handle' ) && strncmp( func2str( obj ), '@(', 2 );
end
