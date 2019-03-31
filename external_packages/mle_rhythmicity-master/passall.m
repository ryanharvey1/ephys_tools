function [ varargout ] = passall( fun, args )
%PASSALL Passes all vector of values to a function
%   
% INPUT
%   fun: A function handle
%   args: Arguments of the function in a vector
%
% Written 2015-11-04 by Jason R. Climer (jason.r.climer@gmail.com)
if ~iscell(args), args = num2cell(args); end
varargout = cell(max(nargout,1),1);
[ varargout{:} ] = fun(args{:});

end

