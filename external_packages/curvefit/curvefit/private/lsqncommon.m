function [x,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB,msg] = lsqncommon(FUN,x,YDATA,LB,UB,options,defaultopt,caller,computeLambda,varargin)
%LSQNCOMMON Solves non-linear least squares problems.
%   [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN]=...
%      LSQNCOMMON(FUN,X0,YDATA,LB,UB,OPTIONS,DEFAULTOPT,CALLER,COMPUTELAMBDA,XDATA,VARARGIN...) 
%   contains all the setup code common to both LSQNONLIN and LSQCURVEFIT to call either the 
%   large-scale SNLS or the medium-scale NLSQ.

%   Copyright 2001-2011 The MathWorks, Inc.
%   $Revision: 1.11.2.15 $  $Date: 2011/05/09 00:39:30 $

% Note: none of these warnings should go to the command line. They've
% all been taken care of in FIT. If they appear from here (not FIT), 
% it is an error.

xstart=x(:);
numberOfVariables=length(xstart);
% Note: XDATA is bundled with varargin already for lsqcurvefit 
lenVarIn = length(varargin);

large = 'large-scale';
medium = 'medium-scale';

% Options setup
switch optimget(options,'Display',defaultopt,'fast')
case {'notify'}
   verbosity = 1;
case {'off','none'}
   verbosity = 0;
case 'iter'
   verbosity = 3;
case 'final'
   verbosity = 2;
case 'testing'
    verbosity = Inf;
otherwise
   verbosity = 2;
end

l = LB; u = UB;
lFinite = ~isinf(l);
uFinite = ~isinf(u);

if min(min(u-xstart),min(xstart-l)) < 0
    xstart = startx(u,l); 
end

gradflag =  strcmp(optimget(options,'Jacobian',defaultopt,'fast'),'on');
line_search = strcmp(optimget(options,'LargeScale',defaultopt,'fast'),'off'); % 0 means large-scale, 1 means medium-scale
mtxmpy = optimget(options,'JacobMult',[]); % use old
if ~isempty(mtxmpy) && ( (ischar(mtxmpy) && isequal(lower(mtxmpy),'atamult')) || ...
                       (isequal(mtxmpy, @atamult)) )
    warning(message('curvefit:lsqncommon:nameClash'));
end

% Convert to inline function as needed
if ~isempty(FUN)  % will detect empty string, empty matrix, empty cell array
    [funfcn] = fprefcnchk(FUN,caller,lenVarIn,gradflag);
else
    invalidFunctionError();
end

fuser = [];  
JAC = [];
x(:) = xstart;
try
    switch funfcn{1}
    case 'fun'
        fuser = feval(funfcn{3},x,varargin{:});
    case 'fungrad'
        [fuser,JAC] = feval(funfcn{3},x,varargin{:});
    case 'fun_then_grad'
        fuser = feval(funfcn{3},x,varargin{:}); 
        JAC = feval(funfcn{4},x,varargin{:});
    otherwise
        error(message('curvefit:lsqncommon:undefinedCalltype', upper( caller )));
    end
catch fevalException
    if ~(isequal(funfcn{1},'fun_then_grad')) || isempty(fuser)
        badfunfcn = funfcn{3};
    else % 'fun_then_grad && ~isempty(fuser) (so it error'ed out on the JAC call)
        badfunfcn = funfcn{4};
    end   
    if isa(badfunfcn,'inline')
        errmsg = sprintf(['User supplied expression or inline function ==> %s\n' ...
                'failed with the following error:\n\n%s'],...
                 formula(badfunfcn),fevalException.message);
        errid = 'curvefit:lsqncommon:userExpOrInlineFcnFailed'; 
    else % function (i.e. string name of), function handle, inline function, or other
        % Convert to char if possible
        try
            charOfFunction = char(badfunfcn);
        catch ignore %#ok<NASGU>
            charOfFunction = '';
        end
        if ~isempty(charOfFunction)
            errmsg = sprintf(['User supplied function ==> %s\n' ...
                    'failed with the following error:\n\n%s'],...
                     charOfFunction,fevalException.message);
        else
            errmsg = sprintf(['User supplied function\n' ...
                    'failed with the following error:\n\n%s'],...
                    fevalException.message);
        end
        errid = 'curvefit:lsqncommon:userFcnFailed'; 
    end
    % Throw an exception with the computed message, but also pass along the
    % original cause of the exception to the caller.
    lsqnException = MException( errid, '%s', errmsg );
    lsqnException = addCause( lsqnException, fevalException );
    throw( lsqnException );
end

if isequal(caller,'lsqcurvefit')
    if ~isequal(size(fuser), size(YDATA))
        error(message('curvefit:lsqncommon:incommensurateValueAndSize'))
    end
    fuser = fuser - YDATA;  % preserve fuser shape until after subtracting YDATA 
end

f = fuser(:);
nfun=length(f);

if gradflag
    % check size of JAC
    [Jrows, Jcols]=size(JAC);
    if isempty(mtxmpy) 
        % Not using 'JacobMult' so Jacobian must be correct size
        if Jrows~=nfun || Jcols ~=numberOfVariables
            error(message('curvefit:lsqncommon:incorrectJacobianSize', nfun, numberOfVariables));
        end
    end
else
    Jrows = nfun; 
    Jcols = numberOfVariables;   
end

% trustregion and enough equations (as many as variables) 
if ~line_search && nfun >= numberOfVariables 
    OUTPUT.algorithm = large;
    
    % trust region and not enough equations -- switch to line_search
elseif ~line_search && nfun < numberOfVariables 
    warning(message('curvefit:lsqncommon:notEnoughVariablesForLargeScale'));
    OUTPUT.algorithm = medium;
    
    % line search and no bounds  
elseif line_search && isempty(l(lFinite)) && isempty(u(uFinite))
    OUTPUT.algorithm = medium;
    
    % line search and  bounds  and enough equations, switch to trust region 
elseif line_search && (~isempty(l(lFinite)) || ~isempty(u(uFinite))) && nfun >= numberOfVariables
    warning(message('curvefit:lsqncommon:invalidBoundsForLMandGN'));
    OUTPUT.algorithm = large;
    
    % can't handle this one:   
elseif line_search && (~isempty(l(lFinite)) || ~isempty(u(uFinite)))  && nfun < numberOfVariables
    error(message('curvefit:lsqncommon:abortingFit'));
end

% Execute algorithm
if isequal(OUTPUT.algorithm,large)
    
    % TODO: The large scale algorithm is done else where. We should not hit this
    % code? Should we warn? error?
    
    if ~gradflag % provide sparsity of Jacobian if not provided.
        Jstr = optimget(options,'JacobPattern',[]);
        if isempty(Jstr)  
            % Put this code separate as it might generate OUT OF MEMORY error
            Jstr = sparse(ones(Jrows,Jcols));
        end
        if ischar(Jstr) 
            if isequal(lower(Jstr),'sparse(ones(jrows,jcols))')
                Jstr = sparse(ones(Jrows,Jcols));
            else
                error(message('curvefit:lsqncommon:JacobPatternMustBeMatrix'));
            end
        end
    else
        Jstr = [];
    end
    warnstate = warning('off', 'all');
    try
        [x,FVAL,LAMBDA,JACOB,EXITFLAG,OUTPUT,msg] = ...
            snls(funfcn,x,l,u,verbosity,options,defaultopt, ...
                 f,JAC,YDATA,caller,Jstr,computeLambda,varargin{:});
    catch e
         warning(warnstate);
         rethrow( e );
    end
    warning(warnstate);
    % So nlsq and snls return the same thing when no bounds exist:
    if all(~isfinite(l)) && all(~isfinite(u))
        LAMBDA.upper=[]; LAMBDA.lower=[];   
    end        
else
    warnstate = warning('off', 'all');
    try
        [x,FVAL,JACOB,EXITFLAG,OUTPUT,msg] = ...
            nlsq(funfcn,x,verbosity,options,defaultopt, ...
                 f,JAC,YDATA,caller,varargin{:});
    catch e
        warning(warnstate);
        rethrow( e );
    end
    warning(warnstate);
    LAMBDA.upper = zeros( size( x ) );
    LAMBDA.lower = zeros( size( x ) );
end
Resnorm = FVAL'*FVAL;
if verbosity > 1 || ( verbosity > 0 && EXITFLAG <=0 )
    disp(msg);
end


% Reset FVAL to shape of the user-function, fuser
FVAL = reshape(FVAL,size(fuser));

%--end of lsqncommon--

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [allfcns,msg] = fprefcnchk(funstr,caller,lenVarIn,gradflag)
%PREFCNCHK Pre- and post-process function expression for FUNCHK.
%   [ALLFCNS,MSG] = PREFUNCHK(FUNSTR,CALLER,lenVarIn,GRADFLAG) takes
%   the (nonempty) expression FUNSTR from CALLER with LenVarIn extra arguments,
%   parses it according to what CALLER is, then returns a string or inline
%   object in ALLFCNS.  If an error occurs, this message is put in MSG.
%
%   ALLFCNS is a cell array: 
%    ALLFCNS{1} contains a flag 
%    that says if the objective and gradients are together in one function 
%    (calltype=='fungrad') or in two functions (calltype='fun_then_grad')
%    or there is no gradient (calltype=='fun'), etc.
%    ALLFCNS{2} contains the string CALLER.
%    ALLFCNS{3}  contains the objective function
%    ALLFCNS{4}  contains the gradient function (transpose of Jacobian).
%  
%    NOTE: we assume FUNSTR is nonempty.
% Initialize
allfcns = {};
gradfcn = [];

if gradflag
    calltype = 'fungrad';
else
    calltype = 'fun';
end

% {fun}
if isa(funstr, 'cell') && length(funstr)==1
    % take the cellarray apart: we know it is nonempty
    if gradflag
        calltype = 'fungrad';
    end
    [funfcn, msg] = fcnchk(funstr{1},lenVarIn);
    if ~isempty(msg)
        error(msg);
    end
    
    % {fun,[]}      
elseif isa(funstr, 'cell') && length(funstr)==2 && isempty(funstr{2})
    if gradflag
        calltype = 'fungrad';
    end
    [funfcn, msg] = fcnchk(funstr{1},lenVarIn);
    if ~isempty(msg)
        error(msg);
    end  
    
    % {fun, grad}   
elseif isa(funstr, 'cell') && length(funstr)==2 % and ~isempty(funstr{2})
    
    [funfcn, msg] = fcnchk(funstr{1},lenVarIn);
    if ~isempty(msg)
        error(msg);
    end  
    [gradfcn, msg] = fcnchk(funstr{2},lenVarIn);
    if ~isempty(msg)
        error(msg);
    end
    calltype = 'fun_then_grad';
    if ~gradflag
        warning(message('curvefit:lsqncommon:invalidFitoptions'));
        calltype = 'fun';
    end   
    
elseif ~isa(funstr, 'cell')  %Not a cell; is a string expression, function name string, function handle, or inline object
    [funfcn, msg] = fcnchk(funstr,lenVarIn);
    if ~isempty(msg)
        error(msg);
    end   
    if gradflag % gradient and function in one file
        gradfcn = funfcn; % Do this so graderr will print the correct name
    end  
else
    invalidFunctionError();
end

allfcns{1} = calltype;
allfcns{2} = caller;
allfcns{3} = funfcn;
allfcns{4} = gradfcn;
allfcns{5}=[];


function invalidFunctionError()
error(message('curvefit:lsqncommon:invalidFUN'))
