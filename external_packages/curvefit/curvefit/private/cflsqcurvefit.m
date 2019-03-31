function [x, Resnorm, FVAL, EXITFLAG, OUTPUT, LAMBDA, JACOB, msg] = cflsqcurvefit( ...
    FUN, x, XDATA, YDATA, LB, UB, fitopt, varargin)
%CFLSQCURVEFIT   Solves non-linear least squares problems.
%
%   [X, RESNORM, FVAL, EXITFLAG, OUTPUT, LAMBDA, JACOB, MSG] = CFLSQCURVEFIT( ...
%       FUN, X, XDATA, YDATA, LB, UB, FITOPT, ...)
%
%   See also FIT, LSQCURVEFIT.

%   Copyright 2001-2010 The MathWorks, Inc.
%   $Revision: 1.10.2.13 $  $Date: 2011/03/28 04:28:38 $

if strcmpi( fitopt.Algorithm, 'gauss-newton' )
    % 'Gauss-Newton' will be removed in a future release. Use
    % 'Levenberg-Marquardt' as the value for 'Algorithm' instead. 
    [x,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB,msg] = iGaussNewton( ...
        FUN,x,XDATA,YDATA,LB,UB,fitopt,varargin{:});
    return
end

% Convert FITOPTIONS to OPTIMOPTIONS
options = optimset( rmfield( getallfields( fitopt ), 'Algorithm' ) );

% Set the options related to the different optimization algorithms
switch lower(fitopt.Algorithm)
    case 'trust-region'
        options.Algorithm = 'trust-region-reflective';
        Jrows = length( XDATA );
        Jcols = length( x );
        options.JacobPattern = ones( Jrows, Jcols );
        options.PrecondBandWidth = Inf;
        iterDisplayFcn = @iIterDispOutFcnTRR;
        
    case 'levenberg-marquardt'
        options.Algorithm = 'levenberg-marquardt';
        iterDisplayFcn = @iIterDispOutFcnLM;
        
    case 'gauss-newton'
        error(message('curvefit:cflsqcurvefit:GaussNewtonIsDeprecated:Error'));
    otherwise
        error(message('curvefit:cflsqcurvefit:invalidAlgorithm'));
end

% Set the interrupt function
options.OutputFcn = @iFittingInterrupt;

% Turn off the display of messages that are built in
options.Display = 'off';
% But turn on our own iterative display if requested
if strcmpi( fitopt.Display, 'iter' )
    options.OutputFcn = {options.OutputFcn, iterDisplayFcn};
end

% Call LSCFTSH
[x,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB] = cflscftsh(FUN,x,XDATA,YDATA,LB,UB,options,varargin{:});

% Generate and display exit message
msg = iExitMessageFromExitFlag( EXITFLAG );
OUTPUT.message = msg;

isDisplayExitMsg = strcmpi( fitopt.Display, 'iter' ) || strcmpi( fitopt.Display, 'final' );
if isDisplayExitMsg 
    disp( msg )
end
end

function stop = iFittingInterrupt(~, ~, ~, varargin)
% iOutputFcn(x,optimValues,state)
stop = cfInterrupt( 'get' );
end

function stop = iIterDispOutFcnTRR(~, optimValues, state, varargin)
%iIterDispOutFcnTrr Output function that produces iterative display.
%
% Helper function that produces iterative display for the trust-region-reflective
% algorithm in lsqcurvefit.

% This output function only displays iterative output; it never stops the run
stop = false;

header = sprintf(['\n                                         Norm of      First-order \n',...
    ' Iteration  Func-count     f(x)          step          optimality   CG-iterations']);
% Format for 0th iteration
formatstrFirstIter = ' %5.0f      %5.0f   %13.6g                  %12.3g';
% Format for iterations >= 1
formatstr = ' %5.0f      %5.0f   %13.6g  %13.6g   %12.3g      %7.0f';

if strcmpi(state,'init')
    % Iterative display header
    disp(header);
elseif strcmpi(state,'iter')
    % Iterative display
    if optimValues.iteration == 0
        % Only a subset of the displayed quantities are available at iteration zero
        fprintf([formatstrFirstIter '\n'], ...
            optimValues.iteration, ...
            optimValues.funccount, ...
            optimValues.resnorm, ...
            norm(optimValues.firstorderopt,Inf));
    else
        fprintf([formatstr '\n'], ...
            optimValues.iteration, ...
            optimValues.funccount, ...
            optimValues.resnorm, ...
            optimValues.stepsize, ...
            norm(optimValues.firstorderopt,Inf), ...
            optimValues.cgiterations);
    end
end
end

function stop = iIterDispOutFcnLM(~, optimValues, state, varargin)
%iIterDispOutFcnLm Output function that produces iterative display.
%
% Helper function that produces iterative display for the levenberg-marquardt
% algorithm in lsqcurvefit.

% This output function only displays iterative output; it never stops the run
stop = false;

if strcmpi(state,'init')
    % Iterative display header
    fprintf( ...
        ['\n                                        First-Order                    Norm of \n', ...
        ' Iteration  Func-count    Residual       optimality      Lambda           step\n']);
elseif strcmpi(state,'iter')
    % Iterative display
    fprintf(' %5.0f       %5.0f   %13.6g    %12.3g %12.6g   %12.6g\n',optimValues.iteration, ...
        optimValues.funccount,optimValues.resnorm,norm(optimValues.gradient,Inf),optimValues.lambda, ...
        norm(optimValues.searchdirection));
end
end

function msg = iExitMessageFromExitFlag( flag )
% Generate an exit message from an exit flag

switch flag
    case 1
        msg = xlate( 'Success. Fitting converged to a solution.' );
    case 2
        msg = xlate( 'Success, but fitting stopped because change in coefficients less than tolerance (TolX).' );
    case 3
        msg = xlate( 'Success, but fitting stopped because change in residuals less than tolerance (TolFun).' );
    case 4
        msg = xlate( 'Success, but fitting stopped because magnitude of search direction smaller than tolerance (TolX).' );
    case 0
        msg = xlate( 'Fitting stopped because the number of iterations or function evaluations exceeded the specified maximum.' );
    case -1
        msg = xlate( 'Fitting stopped by user.' );
    case -2
        msg = xlate( 'Fit not computed because lower bounds are greater than upper bounds.' );
    case -3
        msg = xlate( 'Fitting stopped because the Levenberg-Marquardt regularization parameter became too large.' );
    case -4
        msg = xlate( 'Fitting could not make further progress.' );
    otherwise
        % We should not hit this case.
        msg = xlate( 'Fitting stopped for unknown reason.' );
        warning(message('curvefit:cflsqcurvefit:UnknownExitFlag'));
end
end

function [x,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB,msg] = iGaussNewton( ...
    FUN,x,XDATA,YDATA,LB,UB,fitopt,varargin)
% This is deprecated code to run the Gauss-Newton solver

% ------------Initialization----------------
defaultopt = struct('Display','final','LargeScale','on', ...
   'TolX',1e-6,'TolFun',1e-6,'DerivativeCheck','off',...
   'Diagnostics','off',...
   'Jacobian','off','JacobMult',[],...% JacobMult set to [] by default
   'JacobPattern','sparse(ones(Jrows,Jcols))',...
   'MaxFunEvals','100*numberOfVariables',...
   'DiffMaxChange',1e-1,'DiffMinChange',1e-8,...
   'PrecondBandWidth',0,'TypicalX','ones(numberOfVariables,1)',...
   'MaxPCGIter','max(1,floor(numberOfVariables/2))', ...
   'TolPCG',0.1,'MaxIter',400); 


% If just 'defaults' passed in, return the default options in X
if nargin==1 && nargout <= 1 && isequal(FUN,'defaults')
    x = defaultopt;
    return
end

if nargin < 7, fitopt = [];
    if nargin < 6, UB = [];
        if nargin < 5, LB = [];
            if nargin < 4,
                error(message('curvefit:cflsqcurvefit:wrongNumArgs'));
            end
        end
    end
end
if nargout > 5
    computeLambda = 1;
else
    computeLambda = 0;
end

% Convert FITOPTIONS to OPTIMOPTIONS
pvlist = {};
% We have already established that
%   lower(fitopt.Algorithm) == 'gauss-newton', hence
pvlist(end+1:end+2) = {'LargeScale','off'};

options = optimset(rmfield(getallfields(fitopt),'Algorithm'), pvlist{:});

% XDATA put at the end so it's bundled with varargin in lsqncommon
caller = 'lsqcurvefit';
[x,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB,msg] = ...
    lsqncommon(FUN,x,YDATA,LB,UB,options,defaultopt,caller,...
    computeLambda,XDATA,varargin{:});

FVAL = -FVAL;  % so that FVAL = YDATA - F(X,...);
end

