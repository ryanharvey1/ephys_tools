function [line1,line2] = makedisplay(obj,objectname)

%   Copyright 1999-2009 The MathWorks, Inc.
%   $Revision: 1.8.2.7 $  $Date: 2010/10/08 16:37:25 $

% If no object name is given then use the dependname.
if isempty( objectname )
    objectname = dependnames( obj );
    objectname = objectname {1};
end

line1 = sprintf('%s =', objectname);

if (obj.isEmpty)
    line2 = 'Model:  (empty)';
else
    switch obj.fCategory
    case 'custom'
        if islinear(obj)
            line2a = sprintf('Linear model:\n     ');
        else
            line2a = sprintf('General model:\n     ');
        end
        line2b = fcnstring( obj, objectname );
    case {'spline','interpolant', 'lowess'}
        line2a = sprintf('%s:\n     ',obj.fTypename);
        line2b = sprintf( '  %s', nonParametricFcnString( obj, objectname, argnames( obj ), obj.coeff(1,:) ) );
    case 'library'
        numargs = obj.numArgs;
        args = obj.args;
        if numargs>10
           args = trimargs(args);
           numargs = size(args,1);
        end
        if islinear(obj)
            line2a = sprintf('Linear model %s:\n     ',obj.fTypename);
        else
            line2a = sprintf('General model %s:\n     ',obj.fTypename);
        end
        line2b = fcnstring( obj, objectname, numargs, args );
    otherwise
        error(message('curvefit:fittype:makedisplay:UnknownType'));
    end
    line2 = sprintf('%s%s',line2a,line2b);
end

% ------------- elide some arguments if there are too many

function args = trimargs(args)

% Arguments may include a1,a2,a3,a4,...
%                    or a1,b1,a2,b2,a3,b3,...
%                    or a1,b1,c1,a2,b2,c2,a3,b3,c3,...

% Find a3
i1 = strmatch('a3',args,'exact');
if isempty(i1), return; end

% Find last-numbered a that is a4 or higher
i2 = [];
for j=4:100
   aname = sprintf('a%d',j);
   k = strmatch(aname,args,'exact');
   if isempty(k)
      break
   else
      i2 = k;
   end
end

% If we found both, and they're separated by more than 1, wipe out stuff
if isempty(i2), return; end
if i2<=i1+1, return; end
args(i1,:) = ' ';
if size(args,2)<3
   args(:,3) = ' ';
end
args(i1,1:3) = '...';
args(i1+1:i2-1,:) = [];
