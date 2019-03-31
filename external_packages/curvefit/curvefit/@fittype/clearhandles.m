function obj = clearhandles(obj)
%CLEARHANDLES Method to remove function handles, called by the
%  SAVEOBJ method for this class or a derived class

%   Copyright 2001-2008 The MathWorks, Inc.
%   $Revision: 1.4.2.2 $  $Date: 2008/10/31 05:56:33 $

% Remove function handles that will be useless if loaded later
if isa(obj.expr,'function_handle')
   obj.expr = [];
end
if isa(obj.derexpr,'function_handle')
   obj.derexpr = [];
end
if isa(obj.intexpr,'function_handle')
   obj.intexpr = [];
end
if isa(obj.fStartpt,'function_handle')
   obj.fStartpt = [];
end

