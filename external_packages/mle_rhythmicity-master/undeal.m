function [ out ] = undeal( fun, nout, varargin )
%UNDEAL Summary of this function goes here
%   Detailed explanation goes here
    out = cell(nout,1);
    [out{:}] = deal(fun(varargin{:}));
end

