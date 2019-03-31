function [x,CostFunction,JAC,EXITFLAG,OUTPUT,msg] = nlsq(funfcn,x,verbosity,options,defaultopt,CostFunction,JAC,YDATA,~,varargin)
%NLSQ Solves non-linear least squares problems.
%   NLSQ is the core code for solving problems of the form:
%   min  sum {FUN(X).^2}    where FUN and X may be vectors or matrices.   
%             x
%

%   Copyright 2001-2010 The MathWorks, Inc.
%   $Revision: 1.12.2.9 $  $Date: 2011/03/28 04:28:39 $

%   The default algorithm is the Levenberg-Marquardt method with a 
%   mixed quadratic and cubic line search procedure.  A Gauss-Newton
%   method is selected by setting  OPTIONS.LevenbergMarq='on'. 
%

% ------------Initialization----------------


% Initialization
XOUT = x(:);
% numberOfVariables must be the name of this variable
numberOfVariables = length(XOUT);
msg = '';
how = [];
OUTPUT = [];
EXITFLAG = 1;  %assume convergence
currstepsize = 0;

formatstr=' %5.0f       %5.0f   %13.6g %12.3g %12.3g  ';

% Fixed options (removed from optimset)
lineSearchType = false;  % Use quadcubic line-search
levMarq = false;         % Use Gauss-Newton algorithm

% options
gradflag =  strcmp(optimget(options,'Jacobian',defaultopt,'fast'),'on');
tolX = optimget(options,'TolX',defaultopt,'fast');
tolFun = optimget(options,'TolFun',defaultopt,'fast');
DiffMinChange = optimget(options,'DiffMinChange',defaultopt,'fast');
DiffMaxChange = optimget(options,'DiffMaxChange',defaultopt,'fast');
DerivativeCheck = strcmp(optimget(options,'DerivativeCheck',defaultopt,'fast'),'on');
maxFunEvals = optimget(options,'MaxFunEvals',defaultopt,'fast');
maxIter = optimget(options,'MaxIter',defaultopt,'fast');
if ischar(maxFunEvals)
   if isequal(lower(maxFunEvals),'100*numberofvariables')
      maxFunEvals = 100*numberOfVariables;
   else
      error(message('curvefit:nlsq:maxFunEvalsMustBeInt'))
   end
end


nfun=length(CostFunction);
iter = 0;
numFunEvals = 0;
numGradEvals = 0;
MATX=zeros(3,1);
MATL=[CostFunction'*CostFunction;0;0];
FIRSTF=CostFunction'*CostFunction;
OLDX=lsinit(XOUT,CostFunction,verbosity,levMarq);
PCNT = 0;
EstSum=0.5;
% system of equations or overdetermined
if nfun >= numberOfVariables
   GradFactor = 0;  
else % underdetermined: singularity in JAC'*JAC or GRAD*GRAD' 
   GradFactor = 1;
end
CHG = repmat(1e-7*abs(XOUT)+1e-7, numberOfVariables,1);
status=-1;

while status~=1 && ~(cfInterrupt( 'get' ) && iter >= 1)
   % Interrupt optimization if requested AND done at least one iteration
    
   iter = iter + 1;
   % Work Out Gradients
   if ~(gradflag) || DerivativeCheck
      JACFD = zeros(nfun, numberOfVariables);  % set to correct size
      OLDF=CostFunction;
      CHG = sign(CHG+eps).*min(max(abs(CHG),DiffMinChange),DiffMaxChange);
      for gcnt=1:numberOfVariables
         temp = XOUT(gcnt);
         XOUT(gcnt) = temp +CHG(gcnt);
         x(:) = XOUT;
         CostFunction = feval(funfcn{3},x,varargin{:});  
         if ~isempty(YDATA)
            CostFunction = CostFunction - YDATA;
         end
         CostFunction = CostFunction(:);

         JACFD(:,gcnt)=(CostFunction-OLDF)/(CHG(gcnt));
         XOUT(gcnt) = temp;
      end
      CostFunction = OLDF;
      numFunEvals=numFunEvals+numberOfVariables;
      % Gradient check
      if DerivativeCheck == 1 && gradflag
         
         if isa(funfcn{3},'inline') 
            % if using inlines, the gradient is in funfcn{4}
            graderr(JACFD, JAC, formula(funfcn{4})); %
         else 
            % otherwise fun/grad in funfcn{3}
            graderr(JACFD, JAC,  funfcn{3});
         end
         DerivativeCheck = 0;
      else
         JAC = JACFD;
      end
   else
      x(:) = XOUT;
   end
   
   % Try to set difference to 1e-8 for next iteration
   if nfun==1  % JAC is a column vector or scalar, don't sum over rows
      sumabsJAC = abs(JAC);
   else  % JAC is a row vector or matrix
      sumabsJAC = sum(abs(JAC)')';  % Sum over the rows of JAC
   end
   nonZeroSum = (sumabsJAC~=0);
   CHG(~nonZeroSum) = Inf;  % to avoid divide by zero error
   CHG(nonZeroSum) = nfun*1e-8./sumabsJAC(nonZeroSum);
   
   GradF = 2*(CostFunction'*JAC)'; %2*GRAD*CostFunction;
   NewF = CostFunction'*CostFunction;
   %---------------Initialization of Search Direction------------------
   if status==-1
      if cond(JAC)>1e8 
%         SD=-(GRAD*GRAD'+(norm(GRAD)+1)*(eye(numberOfVariables,numberOfVariables)))\(GRAD*CostFunction);
         SD=-(JAC'*JAC+(norm(JAC)+1)*(eye(numberOfVariables,numberOfVariables)))\(CostFunction'*JAC)';
         if levMarq 
            GradFactor=norm(JAC)+1; 
         end
         how='COND';
      else
         % SD=JAC\(JAC*X-F)-X;
         SD=-(JAC'*JAC+GradFactor*(eye(numberOfVariables,numberOfVariables)))\(CostFunction'*JAC)';
      end
      FIRSTF=NewF;
      GDOLD=GradF'*SD;
      % currstepsize controls the initial starting step-size.
      % If currstepsize has been set externally then it will
      % be non-zero, otherwise set to 1.
      if currstepsize == 0, 
         currstepsize=1; 
      end
      if verbosity > 2
         fprintf(formatstr,iter,numFunEvals,NewF,currstepsize,GDOLD);
      end
      XOUT=XOUT+currstepsize*SD;
      if levMarq
         newf=JAC*SD+CostFunction;
         GradFactor=newf'*newf;
         SD=-(JAC'*JAC+GradFactor*(eye(numberOfVariables,numberOfVariables)))\(CostFunction'*JAC)'; 
      end
      newf=JAC*SD+CostFunction;
      XOUT=XOUT+currstepsize*SD;
      EstSum=newf'*newf;
      status=0;
      if lineSearchType==0; 
         PCNT=1; 
      end
      
   else
      %-------------Direction Update------------------
      gdnew=GradF'*SD;
      if verbosity > 2 
         num=sprintf(formatstr,iter,numFunEvals,NewF,currstepsize,gdnew);
      end
      if gdnew>0 && NewF>FIRSTF
         % Case 1: New function is bigger than last and gradient w.r.t. SD -ve
         % ... interpolate. 
         how='inter';
         [stepsize]=cubici1(NewF,FIRSTF,gdnew,GDOLD,currstepsize);
         currstepsize=0.9*stepsize;
      elseif NewF<FIRSTF
         %  New function less than old fun. and OK for updating 
         %         .... update and calculate new direction. 
         [newstep,fbest] =cubici3(NewF,FIRSTF,gdnew,GDOLD,currstepsize);
         if fbest>NewF,
            fbest=0.9*NewF; 
         end 
         if gdnew<0
            how='incstep';
            if newstep<currstepsize,  
               newstep=(2*currstepsize+1e-4); how=[how,'IF']; 
            end
            currstepsize=abs(newstep);
         else 
            if currstepsize>0.9
               how='int_step';
               currstepsize=min([1,abs(newstep)]);
            end
         end
         % SET DIRECTION.      
         % Gauss-Newton Method    
         temp=1;
         if ~levMarq
            if currstepsize>1e-8 && cond(JAC)<1e8
               SD=JAC\(JAC*XOUT-CostFunction)-XOUT;
               if SD'*GradF>eps,
                  how='ERROR- GN not descent direction';
               end
               temp=0;
            else
               if verbosity > 0
                  warning(message('curvefit:nlsq:poorGradientConditioning'));
               end
               how='CHG2LM';
               levMarq=1;
               currstepsize=abs(currstepsize);               
            end
         end
         
         if (temp)      
            % Levenberg_marquardt Method N.B. EstSum is the estimated sum of squares.
            %                                 GradFactor is the value of lambda.
            % Estimated Residual:
            if EstSum>fbest
               GradFactor=GradFactor/(1+currstepsize);
            else
               GradFactor=GradFactor+(fbest-EstSum)/(currstepsize+eps);
            end
            SD=-(JAC'*JAC+GradFactor*(eye(numberOfVariables,numberOfVariables)))\(CostFunction'*JAC)'; 
            currstepsize=1; 
            estf=JAC*SD+CostFunction;
            EstSum=estf'*estf;
            if verbosity > 2
               num=[num,sprintf('%12.6g ',GradFactor)]; 
            end
         end
         gdnew=GradF'*SD;
         
         OLDX=XOUT;
         % Save Variables
         FIRSTF=NewF;
         GDOLD=gdnew;    
         
         % If quadratic interpolation set PCNT
         if lineSearchType==0, 
            PCNT=1; MATX=zeros(3,1); MATL(1)=NewF; 
         end
      else 
         % Halve Step-length
         how='Red_Step';
         if NewF==FIRSTF,
            msg = sprintf('Optimization terminated successfully: \n Function value converged (no longer changing).');     
            status=1;
            EXITFLAG = 1;
         else
            currstepsize=currstepsize/8;
            if currstepsize<1e-8
               currstepsize=-currstepsize;
            end
         end
      end
      XOUT=OLDX+currstepsize*SD;
      if isinf(verbosity)
         disp([num,'       ',how])
      elseif verbosity > 2
         disp(num)
      end
      
   end %----------End of Direction Update-------------------
   if lineSearchType==0, 
      PCNT=1; MATX=zeros(3,1);  MATL(1)=NewF; 
   end
   % Check Termination 
   if max(abs(SD))< tolX 
       msg = sprintf('Optimization terminated successfully: \n Search direction less than TolX');     
       status=1; EXITFLAG=1;
   elseif (GradF'*SD) < tolFun && ...
           max(abs(GradF)) < 10*(tolFun+tolX)
       msg = sprintf(['Optimization terminated successfully:\n Gradient in the search direction less than TolFun\n',... 
               ' Gradient less than 10*(TolFun+TolX)']);
       status=1; EXITFLAG=1;
       
   elseif numFunEvals > maxFunEvals
             msg = sprintf(['Maximum number of function evaluations exceeded. Increasing\n',...
                            'MaxFunEvals (in fit options) may allow for a better fit, or \n',...
                            'the current equation may not be a good model for the data.']);
       status=1;
       EXITFLAG = 0;
   elseif iter > maxIter
       msg = sprintf(['Maximum number of iterations exceeded. Increasing MaxIter\n',...
                      '(in fit options) may allow for a better fit, or the current\n',...
                      'equation may not be a good model for the data.']);
       status=1;
       EXITFLAG = 0;
   else
       % Line search using mixed polynomial interpolation and extrapolation.
      if PCNT~=0
         % initialize OX and OLDF 
         OX = XOUT; OLDF = CostFunction;
         while PCNT > 0 && ~cfInterrupt( 'get' ) && numFunEvals <= maxFunEvals
            x(:) = XOUT; 
            CostFunction = feval(funfcn{3},x,varargin{:});
            if ~isempty(YDATA)
               CostFunction = CostFunction - YDATA;
            end
            CostFunction = CostFunction(:);
            numFunEvals=numFunEvals+1;
            NewF = CostFunction'*CostFunction;
            % <= used in case when no improvement found.
            if NewF <= OLDF'*OLDF, 
               OX = XOUT; 
               OLDF=CostFunction; 
            end
            [PCNT,MATL,MATX,steplen,NewF,how]=searchq(PCNT,NewF,OLDX,MATL,MATX,SD,GDOLD,currstepsize,how);
            currstepsize=steplen;
            XOUT=OLDX+steplen*SD;
            if NewF==FIRSTF,  
               PCNT=0; 
            end
         end % end while

         XOUT = OX;
         if numFunEvals>maxFunEvals 
             msg = sprintf(['Maximum number of function evaluations exceeded. Increasing\n',...
                            'MaxFunEvals (in fit options) may allow for a better fit, or \n',...
                            'the current equation may not be a good model for the data.']);
            status=1; 
            EXITFLAG = 0;
         end
      end % if PCNT~=0
   end % if max
   
   x(:)=XOUT; 
   switch funfcn{1}
   case 'fun'
      CostFunction = feval(funfcn{3},x,varargin{:});
      if ~isempty(YDATA)
         CostFunction = CostFunction - YDATA;
      end
      CostFunction = CostFunction(:);
      nfun=length(CostFunction);
      % JAC will be updated when it is finite-differenced
   case 'fungrad'
      [CostFunction,JAC] = feval(funfcn{3},x,varargin{:});
      if ~isempty(YDATA)
         CostFunction = CostFunction - YDATA;
      end
      CostFunction = CostFunction(:);
      numGradEvals=numGradEvals+1;
   case 'fun_then_grad'
      CostFunction = feval(funfcn{3},x,varargin{:}); 
      if ~isempty(YDATA)
         CostFunction = CostFunction - YDATA;
      end
      CostFunction = CostFunction(:);
      JAC = feval(funfcn{4},x,varargin{:});
      numGradEvals=numGradEvals+1;
   otherwise
      error(message('curvefit:nlsq:undefinedCalltype'));
   end
   numFunEvals=numFunEvals+1;
      
end  % while
if cfInterrupt( 'get' )
    % Output function terminated the algorithm
    EXITFLAG = -1;
    msg = xlate( 'Fitting stopped by user.' );
end

XOUT=OLDX;
x(:)=XOUT;

OUTPUT.iterations = iter;
OUTPUT.funcCount = numFunEvals;
OUTPUT.firstorderopt = norm(JAC'*CostFunction,inf);

if levMarq
   OUTPUT.algorithm='Levenberg-Marquardt with line search';
else
   OUTPUT.algorithm='Gauss-Newton with line search';
end

%--end of leastsq--


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xold,fold]=lsinit(xnew,fnew,verbosity,levMarq)
%LSINT  Function to initialize NLSQ routine.

xold=xnew;
fold=fnew;



if verbosity>2
   
   if ~levMarq
      if isinf(verbosity)
         header = sprintf(['\n                                                     Directional \n',...
                             ' Iteration  Func-count    Residual     Step-size      derivative   Line-search']);
      else
         header = sprintf(['\n                                                     Directional \n',...
                             ' Iteration  Func-count    Residual     Step-size      derivative ']);
      end
   else
      if isinf(verbosity)
         header = sprintf(['\n                                                     Directional \n',...
                             ' Iteration  Func-count    Residual     Step-size      derivative   Lambda       Line-search']);
      else
         header = sprintf(['\n                                                     Directional \n',...
                             ' Iteration  Func-count    Residual     Step-size      derivative    Lambda']);
      end
   end
   disp(header)
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [pcnt, matl,matx,stepsize,fnew,how]=searchq(pcnt,fnew,~,matl,matx,~,gdold,stepsize,how)
%SEARCHQ Line search routine for FMINU and LEASTSQ functions.
%   Performs line search procedure for unconstrained and least squares
%   optimization. Uses Quadratic Interpolation.
%   When finished pcnt returns 0.

if pcnt==1
% Case 1: Next point less than initial point. 
%     Increase step-length based on last gradient evaluation
    if fnew<matl(1)
% Quadratic Extrapolation using gradient of first point and 
% values of two other points.
        matl(2)=fnew;
        matx(2)=stepsize;
        newstep=-0.5*gdold*stepsize*stepsize/(fnew-gdold*stepsize-matl(1)+eps);
        if newstep<stepsize,how=[how,'QEF ']; newstep=1.2*stepsize; end
        stepsize=1.2*newstep;
        pcnt=2;
    else
% Case 2: New point greater than initial point. Decrease step-length.
        matl(3)=fnew;
        matx(3)=stepsize;
%Interpolate to get stepsize
        stepsize=max([1e-8*stepsize,-gdold*0.5*stepsize^2/(fnew-gdold*stepsize-matl(1)+eps)]);
        how=[how,'r'];
        pcnt=3;
    end
% Case 3: Last run was Case 1 (pcnt=2) and new point less than 
%     both of other 2. Replace. 
elseif pcnt==2  && fnew< matl(2)
    newstep=cubici2(gdold,[matl(1);matl(2);fnew],[matx(1);matx(2);stepsize]);
    if newstep<stepsize,how=[how, 'CEF ']; end
        matl(1)=matl(2);
        matx(1)=matx(2);
        matl(2)=fnew;
        matx(2)=stepsize;
        stepsize=min([newstep,1])+1.5*stepsize;
        stepsize=max([1.2*newstep,1.2*stepsize]);
        how=[how,'i'];
% Case 4: Last run was Case 2: (pcnt=3) and new function still 
%     greater than initial value.
elseif pcnt==3 && fnew>=matl(1)
    matl(2)=fnew;
    matx(2)=stepsize;
    if stepsize<1e-6
        newstep=-stepsize/2;
        % Give up if the step-size gets too small 
        % Stops infinite loops if no improvement is possible.
        if abs(newstep) < (eps * eps), pcnt = 0; end
    else
        newstep=cubici2(gdold,[matl(1);matl(3);fnew],[matx(1);matx(3);stepsize]);
    end
    matx(3)=stepsize;
    if isnan(newstep), stepsize=stepsize/2; else stepsize=newstep; end
    matl(3)=fnew;
    how=[how,'R'];
% Otherwise must have Bracketed Minimum so do quadratic interpolation.
%  ... having just increased step.
elseif pcnt==2 && fnew>=matl(2)
    matx(3)=stepsize;
    matl(3)=fnew;
    [stepsize]=cubici2(gdold,matl,matx);
    pcnt=4;
% ...  having just reduced step.
elseif  pcnt==3  && fnew<matl(1)
    matx(2)=stepsize;
    matl(2)=fnew;
    [stepsize]=cubici2(gdold,matl,matx);
    pcnt=4;
% Have just interpolated - Check to see whether function is any better 
% - if not replace.
elseif pcnt==4 
    pcnt=0;
    stepsize=abs(stepsize);
% If interpolation failed use old point.
    if fnew>matl(2),
        fnew=matl(2);
        how='f';
        stepsize=matx(2);       
    end
end %if pcnt==1
 


function r = cubici1(f2,f1,c2,c1,dx)
%CUBICI1 Cubicly interpolates 2 points and gradients to estimate minimum.
%
%   This function uses cubic interpolation and the values of two 
%   points and their gradients in order to estimate the minimum of a 
%   a function along a line.

if isinf(f2), f2 = 1/eps; end
z = 3*(f1-f2)/dx+c1+c2;
w = real(sqrt(z*z-c1*c2));
r = dx*((z+w-c1)/(c2-c1+2*w));


function step = cubici2(c,f,x)
%CUBICI2 Determine optimizer step from three points and one gradient.
%   STEP = CUBICI2(c,f,x)
%   Finds the cubic p(x) with p(x(1:3)) = f(1:3) and p'(0) = c.
%   Returns the minimizer of p(x) if it is positive.
%   Calls QUADI if the minimizer is negative.

% p(x) = a/3*x^3 - b*x^2 + c*x + d.
% c = p'(0) is the first input parameter.
% Solve [1/3*x.^3 -1/2*x^2 ones(3,1)]*[a b d]' = f - c*x.
% Compute a and b; don't need d.
%    a = 3*(x1^2*(f2-f3) + x2^2*(f3-f1) + x3^2*(f1-f2))/h
%    b = (x1^3*(f2-f3) + x2^3*(f3-f1) + x3^3*(f1-f2))/h
%    where h = (x1-x2)*(x2-x3)*(x3-x1)*(x1*x2 + x2*x3 + x3*x1).
% Local min and max where p'(s) = a*s^2 - 2*b*s + c = 0
% Local min always comes from plus sign in the quadratic formula.
% If p'(x) has no real roots, step = b/a.
% If step < 0, use quadi instead.

x = x(:);
f = f(:);
g = f - c*x;
g = g([2 3 1]) - g([3 1 2]);
y = x([2 3 1]);
h = prod(x-y)*(x'*y);
a = 3*(x.^2)'*g/h;
b = (x.^3)'*g/h;

% Find minimizer.
step = (b + real(sqrt(b^2-a*c)))/a;

% Is step acceptable?
if step < 0 || ~isfinite(step)
   step = abs(quadi(x,f));
end
if isnan(step)
   step = x(2)/2;
end



function [s,f] = cubici3(f2,f1,c2,c1,dx)
%CUBICI3  Cubicly interpolates 2 points and gradients to find step and min.
%
%   This function uses cubic interpolation and the values of 
%   two points and their gradients in order to estimate the minimum s of a 
%   a function along a line and returns s and f=F(s);
%

%  The equation is F(s) = a/3*s^3 + b*s^2 + c1*s + f1
%      and F'(s) = a*s^2+2*b*s + c1
%  where we know that 
%          F(0) = f1
%          F'(0) = c1  
%          F(dx) = f2   implies: a/3*dx^3 + b*dx^2 + c1*dx + f1 = f2
%          F'(dx) = c2  implies: a*dx^2+2*b*dx + c1 = c2

if isinf(f2), 
    f2 = 1/eps; 
end
a = (6*(f1-f2)+3*(c1+c2)*dx)/dx^3;
b = (3*(f2-f1)-(2*c1+c2)*dx)/dx^2;
disc = b^2 - a*c1;
if a==0 && b==0 
    % F(s) is linear: F'(s) = c1, which is never zero;
    % minimum is s=Inf or -Inf (we return s always positive so s=Inf).
    s = inf; 
elseif a == 0 
    % F(s) is quadratic so we know minimum s
    s = -c1/(2*b);
elseif disc <= 0
    % If disc = 0 this is exact. 
    % If disc < 0 we ignore the complex component of the root.
    s = -b/a;  
else
    s = (-b+sqrt(disc))/a;
end
if s<0,  s = -s; end
if isinf(s)
    f = inf;
else
    % User Horner's rule
    f = ((a/3*s + b)*s + c1)*s + f1;
end



function step = quadi(x,f)
%QUADI Determine optimizer step from three points.
%   STEP = QUADI(x,f)
%   Finds the quadratic p(x) with p(x(1:3)) = f(1:3).
%   Returns the minimizer (or maximizer) of p(x).

% p(x) = a*x^2 + b*x + c.
% Solve [x.^2 x ones(3,1)]*[a b c]' = f.
% Minimum at p'(s) = 0,
% s = -b/(2*a) = (x1^2*(f2-f3) + x2^2*(f3-f1) + x3^2*(f1-f2))/ ...
%                 (x1*(f2-f3) + x2*(f3-f1) + x3*(f1-f2))/2

x = x(:); 
f = f(:);
g = f([2 3 1]) - f([3 1 2]);
step = ((x.*x)'*g)/(x'*g)/2;
