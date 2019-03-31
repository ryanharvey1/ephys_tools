function r = iqrCodeGen(x) %#codegen
%IQRCODEGEN Estimate interquartile range 
%   iqrCodeGen returns the interquartile range of the data x, a single-
%   or double-precision vector.
r = iqr(x);
end