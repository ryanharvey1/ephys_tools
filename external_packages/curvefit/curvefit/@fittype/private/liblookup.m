function obj = liblookup(libname,obj)
% LIBLOOKUP set up a library function.
% OBJ_OUT = LIBLOOKUP(FUNNAME, OBJ) initializes OBJ to be the
% library function FUNNAME and returns the result OBJ_OUT.
%
% Note: .expr, .derexpr, .intexpr functions may assume X is a column vector.


%   Copyright 1999-2009 The MathWorks, Inc.
%   $Revision: 1.30.2.8 $  $Date: 2010/10/08 16:37:24 $

switch libname(1:3)
case 'exp'  % exp1 and exp2
    obj.fFeval = 1;
    n = str2double(libname(end))+length(libname)-4;
    obj.coeff = char((0:2*n-1)'+'a');
    obj.numCoeffs = size(obj.coeff,1);
    if n == 1
        obj.defn = 'a*exp(b*x)';
        obj.fNonlinearcoeffs = 2; % b
    elseif n == 2
        obj.defn = 'a*exp(b*x) + c*exp(d*x)';
        obj.fNonlinearcoeffs = [2,4]; % b and d
    else
        error(message('curvefit:fittype:liblookup:nameNotFound', libname));
    end
    obj.fFitoptions = fitoptions('method','nonlinearleastsquares','Jacobian','on');
case 'pow' % power1,power2
    obj.fFeval = 1;
    n = str2double(libname(end))+length(libname)-6;
    if n == 1
        obj.coeff = strvcat('a','b');
        obj.defn = 'a*x^b';
        obj.fNonlinearcoeffs = 2; % b
    elseif n == 2;
        obj.coeff = strvcat('a','b','c');
        obj.defn = 'a*x^b+c';
        obj.fNonlinearcoeffs = 2; % b
    end
    obj.numCoeffs = size(obj.coeff,1);
    obj.fFitoptions = fitoptions('method','nonlinearleastsquares','Jacobian','on');
case 'gau'
    obj.fFeval = 1;
    n = str2double(libname(end))+length(libname)-6;
    if n < 1 || n > 8
        error(message('curvefit:fittype:liblookup:nameNotFound', libname));
    end
    obj.coeff = [repmat(['a';'b';'c';],n,1) ...
            reshape(repmat(num2str((1:n)'),1,3)',3*n,1)];
    obj.numCoeffs = size(obj.coeff,1);
    obj.fNonlinearcoeffs = reshape([2:3:3*n; 3:3:3*n],1,2*n);
    if n > 2
        pstr = sprintf('\n              ');
    else
        pstr = ' ';
    end
    for i = 1:n
        si = num2str(i);
        pstr = sprintf('%sa%s*exp(-((x-b%s)/c%s)^2) + ',pstr,si,si,si);
        if mod(i,2) == 0 && i ~= n
            pstr = sprintf('%s\n              ',pstr);
        end
    end
    pstr(end-2:end)=[];
    obj.defn = pstr;
    obj.fConstants = {n};
    lowerbnds = repmat([-inf;-inf;0],n,1);
    obj.fFitoptions = fitoptions('method','nonlinearleastsquares','Jacobian','on','lower',lowerbnds);
case 'sin'
    obj.fFeval = 1;
    n = str2double(libname(end))+length(libname)-4;
    if n < 1 || n > 9
        error(message('curvefit:fittype:liblookup:nameNotFound', libname));
    end
    obj.coeff = [repmat(['a';'b';'c';],n,1) ...
            reshape(repmat(num2str((1:n)'),1,3)',3*n,1)];
    obj.numCoeffs = size(obj.coeff,1);
    % obj.fNonlinearcoeffs = reshape([2:3:3*n; 3:3:3*n],1,2*n);  
    obj.fNonlinearcoeffs = []; % Seems to work better as not separable
    if n > 3 
        k = 2; 
    else 
        k = 2-n; 
    end
    pstr = ' ';
    for i = 1:n
        if mod(k+i,3) == 0
            pstr = sprintf('%s\n                    ',pstr);
        end
        si = num2str(i);
        pstr = sprintf('%sa%s*sin(b%s*x+c%s) + ',pstr,si,si,si);
    end
    obj.defn = pstr(1:end-3);
    obj.fConstants = {n};
    lowerbnds = repmat([-inf;0;-inf],n,1);
    upperbnds = repmat([inf; inf; inf],n,1);
    obj.fFitoptions = fitoptions('method','nonlinearleastsquares','Jacobian','on',...
        'lower',lowerbnds,'upper',upperbnds);
case 'rat'
    m = str2num(libname(end));
    n = str2num(libname(end-1));
    if m > 5 || n > 5 || length(libname) ~= 5
        error(message('curvefit:fittype:liblookup:nameNotFound', libname));
    end
    obj.fFeval = 1;
    obj.coeff = [repmat('p',n+1,1) num2str((1:n+1)');...
            repmat('q',m,1) num2str((1:m)')]; % n&m < 10;
    obj.numCoeffs = size(obj.coeff,1);
    % obj.fNonlinearcoeffs = [n+2:n+2+m-1];
    obj.fNonlinearcoeffs = []; % Seems to work better as not separable
    if n == 5
        pstr = sprintf('\n               (');
    else
        pstr = '(';
    end
    if (n == 0)
        pstr = sprintf('%sp1)',pstr);
    else % n > 1
        for i = 1 : n-1
            pstr = sprintf('%sp%s*x^%s + ',pstr,num2str(i),num2str(n-i+1));
        end
        pstr = sprintf('%sp%s*x + p%s)',pstr,num2str(n),num2str(n+1));
    end
    if (m == 1)
        qstr = sprintf('(x + q1)');
    else
        qstr = sprintf('(x^%s + ',num2str(m));
        for i = 2 : m-1
            qstr = sprintf('%sq%s*x^%s + ',qstr,num2str(i-1),num2str(m-i+1));
        end
        qstr = sprintf('%sq%s*x + q%s)',qstr,num2str(m-1),num2str(m));
    end
    if (m+n > 4)
        nl = sprintf(' /\n               ');
    else
        nl = sprintf(' / ');
    end
    obj.defn = sprintf('%s%s%s', pstr, nl, qstr);
    obj.fConstants = {n,m};
    obj.fFitoptions = fitoptions('method','nonlinearleastsquares','Jacobian','on');
case 'wei' % weibull a*b*x^(b-1)*e^(-a*x^b)
    obj.fFeval = 1;
    obj.coeff = ['a';'b'];
    obj.numCoeffs = size(obj.coeff,1);
    obj.defn = 'a*b*x^(b-1)*exp(-a*x^b)';
    lowerbnds = [0; 0];
    obj.fFitoptions = fitoptions('method','nonlinearleastsquares','Jacobian','on','lower',lowerbnds);
case 'pol'
    digits = isstrprop( libname, 'digit' );
    if digits(end) && ~digits(end-1)
        obj = iPolynomialCurveLookup( libname, obj );
    elseif digits(end) && digits(end-1)
        obj = iPolynomialSurfaceLookup( libname, obj );
    else
        error(message('curvefit:fittype:liblookup:nameNotFound', libname));
    end
case 'fou'  % Fourier terms
    % Y = A0 + A(1)*cos(x*w)   + B(1)*sin(x*w)+...
    %        ...
    %        + A(n)*cos(n*x*w) + B(n)*sin(n*x*w)
    % where w is the nonlinear parameter to be fitted. Other
    % parameters will be fitted linearly using backslash.

    obj.fFeval = 1;
    n = str2num(libname(end))+length(libname)-8;
    if n < 1 || n > 9
        error(message('curvefit:fittype:liblookup:nameNotFound', libname));
    end
    pstr = strvcat('a0','a1','b1');
    for i = 2:n
        numstr = num2str(i);
        pstr = strvcat(pstr, sprintf('a%s', numstr), sprintf('b%s', numstr));
    end
    pstr = strvcat(pstr, 'w ');
    obj.coeff = pstr; 
    obj.numCoeffs = size(obj.coeff,1);
    obj.fNonlinearcoeffs = obj.numCoeffs;
    if n > 5, 
        pstr = sprintf('\n              ');
    else
        pstr = '';
    end
    pstr = sprintf('%s a0 + a1*cos(x*w) + b1*sin(x*w) + ',pstr);
    if n > 1
      pstr = sprintf('%s\n               ', pstr);
    end
    for i = 2:n
        si = num2str(i);
        pstr = sprintf('%sa%s*cos(%s*x*w) + b%s*sin(%s*x*w) + ',pstr,si,si,si,si);
        if mod(i,2)==1 && i ~= n
            pstr = sprintf('%s\n               ', pstr);
        end
    end
    pstr(end-2:end) = [];

    obj.defn = pstr;
    obj.fConstants = {n};
    obj.fFitoptions = fitoptions('method','nonlinearleastsquares','Jacobian','on');
case {'smo','cub','nea','spl','lin','pch', 'bih'}
    obj.fFeval = 1;
    obj.coeff = 'p';
    obj.numCoeffs = [];
    obj.fFitoptions = iInterpolantFitOptions( libname );
    if size( obj.indep, 1 ) == 1 % curve
        obj.defn = 'piecewise polynomial';
    else % assume size( obj.indep, 1 ) == 2, i.e., a surface
        obj.defn = iSurfaceInterpolantDefinition( libname );
    end
case {'low', 'loe'}
    obj.fFeval = 1;
    obj.coeff = 'p';
    obj.numCoeffs = [];
    if strcmpi( libname, 'lowess' )
        obj.defn = 'lowess (linear) smoothing regression';
    elseif strcmpi( libname, 'loess' )
        obj.defn = 'loess (quadratic) smoothing regression';
    else
        error(message('curvefit:fittype:liblookup:nameNotFound', libname));
    end
    obj.fFitoptions = fitoptions( 'method', 'lowessfit' );
otherwise
    error(message('curvefit:fittype:liblookup:nameNotFound', libname));
end

% Now define function handle fields
obj = sethandles(libname,obj);

end

function opts = iInterpolantFitOptions( libname )
switch libname(1:3)
    case 'smo'
        opts = fitoptions('method','smoothing');
    case 'cub'
        opts = fitoptions('method','cubicsplinei');
    case 'nea'
        opts = fitoptions('method','nearesti');
    case 'spl'
        opts = fitoptions('method','cubicsplinei');
    case 'lin'
        opts = fitoptions('method','lineari');
    case 'pch'
        opts = fitoptions('method','pchipi');
    case 'bih'
        opts = fitoptions('method','BiharmonicInterpolant');
end
end

function defn = iSurfaceInterpolantDefinition( libname )
switch libname(1:3)
    case 'nea' % 'nearestinterp'
        defn = 'piecewise constant surface';
    case 'lin' % 'linearinterp'
        defn = 'piecewise linear surface';
    case 'cub' % 'cubicinterp'
        defn = 'piecewise cubic surface';
    case 'bih' % 'biharmonicinterp '
        defn = 'biharmonic surface';
end
end

function obj = iPolynomialCurveLookup( libname, obj )
% Y = P(1)*X^N + P(2)*X^(N-1) + ... + P(N)*X + P(N+1)
obj.linear = 1;
obj.fFeval = 1;
N = str2num(libname(end))+length(libname)-5;
if N < 1 || N > 9
    error(message('curvefit:fittype:liblookup:nameNotFound', libname));
end
pstr = 'p1 ';
for i = 2:N+1
    numstr = num2str(i);
    if i <= 9
        pstr(i,:) = sprintf('p%s ',numstr);
    else
        pstr(i,:) = sprintf('p%s',numstr);
    end
end
obj.coeff = pstr;  % this should actually be p1,p2, etc
obj.numCoeffs = N+1;
k = floor(N/4);
pstr = '';
for i = 1:N-1
    pstr = sprintf('%sp%s*x^%s + ',pstr,num2str(i),num2str(N-i+1));
    if i+k==6
        pstr = sprintf('%s\n                    ',pstr);
    end
end
obj.defn = sprintf('%sp%s*x + p%s',pstr, num2str(N), num2str(N+1));
obj.fFitoptions = fitoptions('method','linearleastsquares');
end

function obj = iPolynomialSurfaceLookup( libname, obj )
% e.g., p1 + px*x + pxx*x^2 + py*y + pxy*x*y + pxxy*x^2*y
obj.linear = 1;
obj.fFeval = 1;

m = str2double( libname(end-1) );
n = str2double( libname(end)   );
obj.fConstants = {m, n};

maxDegree = max( m, n );
coeffs = {'p00'};
defn = 'p00';
for i = 1:maxDegree
    for j = i:-1:0
        if j <= m && (i-j) <=n
            if j == 0
                tx = '';
            elseif j == 1
                tx = '*x';
            else
                tx = sprintf( '*x^%d', j );
            end
            if (i-j) == 0
                ty = '';
            elseif (i-j) == 1
                ty = '*y';
            else
                ty = sprintf( '*y^%d', (i-j) );
            end
            coeffs{end+1} = sprintf( 'p%d%d', j, (i-j) );
            defn = sprintf( '%s + %s%s%s', defn, coeffs{end}, tx, ty );
        end
    end
end
obj.coeff = strvcat( coeffs{:} );
obj.numCoeffs = length( coeffs );
obj.defn = defn;
obj.fFitoptions = fitoptions( 'method', 'linearleastsquares' );
end
