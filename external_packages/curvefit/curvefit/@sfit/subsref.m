function out = subsref(obj, subs)
%SUBSREF Evaluate SFIT object.

%   Copyright 2008-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2010/10/08 16:37:33 $

if (isempty(obj))
    error(message('curvefit:sfit:subsref:emptyFit'));
end

% In case nested subsref like f.p.coefs
currsubs = subs(1);

switch currsubs.type
    case '()'
        out = iParenthesesReference( obj, currsubs );

    case '.'
        out = iDotReference( obj, currsubs );
        
    otherwise % case '{}'
        error(message('curvefit:sfit:subsref:cellarrayBracketsNotAllowed'))
end

if length(subs) > 1
    subs(1) = [];
    out = subsref(out, subs);
end

end

function out = iParenthesesReference( obj, subs )
inputs = subs.subs;

if (isempty(fevalexpr(obj)))
    out = [];
else
    try
        out = feval( obj, inputs{:} );
    catch ME
        throwAsCaller( ME );
    end
end
end

function out = iDotReference( obj, subs )

argname = subs.subs;
% is it coeff or prob parameter?
coeff = strcmp( argname, coeffnames(obj) );
prob = strcmp( argname, probnames(obj) );

% which index is it?
if any( coeff )
    out = obj.fCoeffValues{coeff};
elseif any( prob )
    out = obj.fProbValues{prob};
else
    % As coefficients and problem parameters must be different, it must be
    % that the name the user gave us is neither coefficient not problem
    % parameter.
    error(message('curvefit:sfit:subsref:invalidName'));
end
end


