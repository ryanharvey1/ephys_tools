function [line1,line2,line3,line4] = makedisplay(obj,objectname,out,clev)
% MAKEDISPLAY

%   Copyright 2008-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2010/10/08 16:37:37 $

if nargin<4, clev = 0.95; end

line1 = sprintf('%s =', objectname);
line4 = ''; %default

if (isempty(obj))
    line2 = xlate('Model: (empty)');
    line3 = xlate('Coefficients: (empty)');
else
    switch category(obj)
        case 'custom'
            if islinear(obj)
                line2a = sprintf('Linear model:\n     ');
            else
                line2a = sprintf('General model:\n     ');
            end
            line2b = fcnstring( obj, objectname, 2, indepnames( obj ) );
            line2c = iMakeScalingLine( obj );
            try
                ci = confint(obj,clev);
                line3a = sprintf('Coefficients (with %g%% confidence bounds):\n',100*clev);
                line3b = argstring(char(coeffnames(obj)),obj.fCoeffValues,...
                    ci,obj.activebounds);
            catch ignore
                line3a = sprintf('Coefficients:\n');
                line3b = argstring(char(coeffnames(obj)),obj.fCoeffValues);
            end
            probnamesarray = char(probnames(obj));
            if ~isempty(probnamesarray)
                line4a = sprintf('Problem parameters:\n');
                line4b = argstring(probnamesarray,obj.fProbValues);
                line4 = sprintf('%s%s',line4a,line4b);
            end
        case {'spline','interpolant', 'lowess'}
            line2a = sprintf('%s:\n     ',prettyname(obj));
            line2b = sprintf( '  %s', nonParametricFcnString( obj, objectname, indepnames( obj ), char(coeffnames(obj)) ) );
            line2c = iMakeScalingLine( obj );
            if nargin>=3 && isfield(out,'p')
                line3a = sprintf('Smoothing parameter:\n');
                line3b = sprintf('       p = %0.8g',out.p);
            else
                line3a = sprintf('Coefficients:\n');
                line3b = argstring(char(coeffnames(obj)),{xlate('coefficient structure')});
            end
        case 'library'
            if islinear(obj)
                line2a = sprintf('Linear model %s:\n     ',prettyname(obj));
            else
                line2a = sprintf('General model %s:\n     ',prettyname(obj));
            end
            line2b = fcnstring(obj, objectname, 2, indepnames( obj ) );
            line2c = iMakeScalingLine( obj );
            try
                ci = confint(obj,clev);
                line3a = sprintf('Coefficients (with %g%% confidence bounds):\n',100*clev);
                line3b = argstring(char(coeffnames(obj)),obj.fCoeffValues,...
                    ci,obj.activebounds);
            catch ignore
                line3a = sprintf('Coefficients:\n');
                line3b = argstring(char(coeffnames(obj)),obj.fCoeffValues);
            end
            probnamesarray = char(probnames(obj));
            if ~isempty(probnamesarray)
                line4a = sprintf('Problem parameters:\n');
                line4b = argstring(probnamesarray,obj.fProbValues);
                line4 = sprintf('%s%s',line4a,line4b);
            end
            
        otherwise
            error(message('curvefit:sfit:makedisplay:unknownFittype'))
    end
    line2 = sprintf('%s%s%s',line2a,line2b,line2c);
    line3 = sprintf('%s%s',line3a,line3b);
end

end

function line = iMakeScalingLine( obj )
indep = indepnames(obj);

xIsScaled = ~isequal( obj.meanx, 0 ) || ~isequal( obj.stdx, 1 );
yIsScaled = ~isequal( obj.meany, 0 ) || ~isequal( obj.stdy, 1 );

if xIsScaled && yIsScaled
    line = sprintf( [
        '\n       where %s is normalized by mean %0.4g and std %0.4g', ...
        '\n       and where %s is normalized by mean %0.4g and std %0.4g'], ...
        indep{1}, obj.meanx, obj.stdx, indep{2}, obj.meany, obj.stdy );

elseif xIsScaled
    line = sprintf('\n       where %s is normalized by mean %0.4g and std %0.4g', ...
        indep{1}, obj.meanx, obj.stdx);
    
elseif yIsScaled
    line = sprintf('\n       where %s is normalized by mean %0.4g and std %0.4g', ...
        indep{2}, obj.meany, obj.stdy );
else 
    line = '';
end

end

