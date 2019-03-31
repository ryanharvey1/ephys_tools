function [n] = plx_close(filename)
import CMBHOME.PLX2Mat.*
% plx_close(filename): Close the .plx file
%
% [n] = plx_close(filename)
%
% INPUT:
%   filename - if empty string, will close any open files
% OUTPUT:
%   n - always 0

[n] = mexPlex(22,filename);
