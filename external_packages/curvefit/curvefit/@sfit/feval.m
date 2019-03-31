function varargout = feval(varargin)
%FEVAL  FEVAL an SFIT object.

%   Copyright 2008 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2008/10/31 05:57:00 $

obj = varargin{1};
if ~isa(obj,'sfit')
    % If any of the elements in varargin are SFIT objects, then the
    %  overloaded SFIT feval is called even if the first argument
    %  is a string.  In this case, we call the builtin feval.
    [varargout{1:max(1,nargout)}] = builtin('feval',varargin{:});
    return
end

% Parse inputs for data
[xdata, ydata] = iXYDataFromInputs( varargin(2:end) );

% Normalize the data
xdata = (xdata - obj.meanx)/obj.stdx;
ydata = (ydata - obj.meany)/obj.stdy;

try
    [varargout{1:max(1,nargout)}] = evaluate(obj, obj.fCoeffValues{:},...
        obj.fProbValues{:}, xdata, ydata);
catch cause
    exception = MException( 'curvefit:sfit:feval:evaluationError', ...
        'Error while trying to evaluate SFIT model.' );
    exception = addCause( exception, cause );
    throw( exception );
end
end

function [xdata, ydata] = iXYDataFromInputs( inputs )

ninputs = length( inputs );
if ninputs < 1
    % FO( X ) --> error!
    throwAsCaller( MException( 'curvefit:sfit:feval:notEnoughInputs', ...
          'Not enough input arguments.\nTo evaluate an SFIT object you must specify either two arguments, or one argument with two columns.') );
      
elseif ninputs == 1
    % FO( XY ) --> X = XY(:,1), Y = XY(:,2)
    if size( inputs{1}, 2 ) == 2
        xdata = inputs{1}(:,1);
        ydata = inputs{1}(:,2);
    else
        throwAsCaller( MException( 'curvefit:sfit:subsref:invalidInput',...
            'To evaluate an SFIT object you must specify either two arguments, or one argument with two columns.') );
    end
    
elseif ninputs == 2
    % FO( X, Y )
    xdata = inputs{1};
    ydata = inputs{2};
    
else % ninputs > 2
    % FO( X, Y, ... ) --> error
    throwAsCaller( MException( 'curvefit:sfit:feval:tooManyInputs', ...
          'Too many input arguments.\nTo evaluate an SFIT object you must specify either two arguments, or one argument with two columns.') );
end

end




