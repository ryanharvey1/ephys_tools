function definition = definitionFromLibname(numindep, libname)
%DEFINITIONFROMLIBNAME   Get definition line from library name
%
%   DEFINITION = DEFINITIONFROMLIBNAME(NUMINDEP, LIBNAME) is a definition (aka
%   formula) suitable for use in the display of a FITTYPE object. DEFINITION
%   contains formats for use with SPRINTF:
%      '%1$s' == left hand side of a defition, e.g., 'cf( x )', 'sf( x, y )'
%      '%2$s' == name of coefficient.
%
%   The library name, LIBNAME, is given by the "fType" property of a FITTYPE
%   object.

%   Copyright 2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $    $Date: 2010/10/08 16:37:23 $ 

switch numindep
    case 1 % curves
        definition = xlate( '%1$s = piecewise polynomial computed from %2$s' );
    case 2 % surfaces
        switch libname
            case 'nearestinterp'
                definition = xlate( '%1$s = piecewise constant surface computed from %2$s' );
            case 'linearinterp'
                definition = xlate( '%1$s = piecewise linear surface computed from %2$s' );
            case 'cubicinterp'
                definition = xlate( '%1$s = piecewise cubic surface computed from %2$s' );
            case 'biharmonicinterp'
                definition = xlate( '%1$s = biharmonic surface computed from %2$s' );
            case 'lowess'
                definition = xlate( '%1$s = lowess (linear) smoothing regression computed from %2$s' );
            case 'loess'
                definition = xlate( '%1$s = loess (quadratic) smoothing regression computed from %2$s' );
            otherwise
                % We shouldn't hit this case
                error(message('curvefit:fittype:InvalidState', libname));
        end
    otherwise
        % We shouldn't hit this case
        error(message('curvefit:fittype:InvalidState', libname));
end

