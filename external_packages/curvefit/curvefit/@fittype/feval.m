function varargout = feval(varargin)
%FEVAL  Evaluate a FITTYPE object.
%   F = FEVAL(FITOBJ,A,B,...,X) evaluates the function value F of FITOBJ with
%   coefficients A,B,... and data X.
%
%   [F, J] = FEVAL(FITOBJ,A,B,...,X) evaluates the function value F and Jacobian
%   J, with respect to the coefficients, of FITOBJ with coefficients A,B,... and
%   data X.
%
%   FEVAL(FITOBJ,A,B,...,X,Y) evaluates FITOBJ with coefficients A,B,... at the
%   data X,Y, where FITOBJ represents a surface (i.e., a function of two
%   variables). 

%   Copyright 1999-2008 The MathWorks, Inc.
%   $Revision: 1.16.2.6 $  $Date: 2008/10/31 05:56:37 $

[varargout{1:nargout}] = evaluate( varargin{:} );
end
