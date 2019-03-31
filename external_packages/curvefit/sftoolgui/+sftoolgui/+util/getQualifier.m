function qualifier = getQualifier(fittype)
%getQualifier returns the "qualifier" of a fittype. 
%   QUALIFIER = getQualifier(FITTYPE) returns the QUALIFIER for the given
%   FITTYPE.
%
%   Examples of qualifiers are the number of polynomial degrees and the
%   number of fourier terms.

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2011/03/22 18:30:11 $

% Get the "type" of the fittype.
theType = type( fittype );

if ~isempty( strfind( theType, 'interp' ) )
    % Interpolants
    qualifier = iGetInterpolantQualifier(theType);
elseif strcmpi('lowess', theType)
    % Lowess
    qualifier = 'linear';
elseif strcmpi('loess', theType)
    % Loess
    qualifier = 'quadratic';
elseif strncmp( 'poly', theType, 4 )
    % Polynomial
    qualifier = iGetPolynomialQualifier(fittype);
elseif strncmpi( theType, 'fourier', 7)
    % Fourier
    qualifier = theType(8);
elseif strncmpi( theType, 'exp', 3)
    % Exponential
    qualifier = theType(4);
elseif strncmpi( theType, 'gauss', 5)
    % Gaussian
    qualifier = theType(6);
elseif strncmpi( theType, 'power', 5)
    % Power
    qualifier = theType(6);
elseif strncmpi( theType, 'rat', 3)
    % Rational
    qualifier = iGetRationalQualifier(fittype);
elseif strncmpi( theType, 'sin', 3)
    % Sum of sine
    qualifier = theType(4);
else
    qualifier = '';
end
end

function qualifier = iGetInterpolantQualifier(theType)

switch theType
    case {'splineinterp', 'cubicinterp'}
        qualifier = 'cubic';
    case 'nearestinterp'
        qualifier = 'nearest';
    case 'linearinterp'
        qualifier = 'linear';
    case 'pchipinterp'
        qualifier = 'PCHIP';
    case 'biharmonicinterp'
        qualifier = 'v4';
    otherwise
        qualifier = 'linear';
end
end

function qualifier = iGetPolynomialQualifier(fittype)
degree = constants( fittype );
if ~isempty(degree)
    qualifier = sprintf('%d%d', degree{:});
else % assume we are dealing with curve models
    theType = type( fittype ); 
    qualifier = theType(5);
end
end

function qualifier = iGetRationalQualifier(fittype)
degree = constants( fittype );
qualifier = sprintf('%d%d', degree{:});
end


