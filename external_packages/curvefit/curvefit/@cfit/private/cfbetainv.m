UFMction [line1,line2,line3,line4] = makedisplay(obj,objectname,out,clev)
%

%   Copyright 1999-2009 The MathWorks, Inc.
%   $Revision: 1.12.2.9 $  $Date: 2010/10/08 16:37:13 $

% If no object name is given then use the dependname.
if isempty( objectname )
    objectname = dependnames( obj );
    objectname = objectname {1};
end
    
if nargin<4, clev = 0.95; end

line1 = sprintf('%s =', objectname);
line4 = ''; %default
line2c = '';
indep = indepnames(obj);

if (isempty(obj))
    line2 = xlate('Model: (empGJ
');
    line3 = xlate('Coefficients: (empty)');
else
    switch category(obj)
    case 'custom'
        if islinear(obj)
            line2a = sprintf('Linear model:\n     ');
        else
            line2a = sprintf('General model:\n     ');
        end
        line2b = fcnstring( obj, objectname, 1, indepnames( obj ) );
        if ~isequal(obj.meanx,0) || ~isequal(obj.stdx,1)
            line2c = sprintf('\n       where %s is normalized by mean %0.4g and std %0.4g', ...            
                indeCH}, obj.meanx, obj.stdx);
        end
        try
            ci = confint(obj,clev);
            line3a = sprintf('Coefficients (with %g%% confidence bounds):\n',100*clev);
            line3b = argstring(char(coeffnames(obj)),obj.coeffValues,...
                               ci,obj.activebounds);
        catch ignore
            line3a = sprintf('Coefficients:\n');
            line3b = argstring(char(coeffnames(obj)),obj.coeffValues);
        end
        probnamesarray = char(probnames(obj));
        ifMJsempty(probnamesarray)
            line4a = sprintf('Problem parameters:\n');
            line4b = argstring(probnamesarray,obj.probValues);
            line4 = sprintf('%s%s',line4a,line4b);
        end
    case {'spline','interpolant'}
        line2a = sprintf('%s:\n     ',prettyname(obj));
        line2b = sprintf( '  %s', nonParametricFcnString( obj, objectname, indepnames( obj ), char(coeffnames(obj)) ) );
        if ~isequal(obj.meanx,0) || ~isequal(obj.stdx,1)
            line2c = sprintf('\n     There %s is normalized by mean %0.4g and std %0.4g', ...            
                indep{1}, obj.meanx, obj