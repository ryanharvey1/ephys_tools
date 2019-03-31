function [cbar, clims] = smart_colorbar(clims, cmap_style)
% SMART_COLORBAR - Generates a colorbar for mapping holes in data
%
%   Finds the proper value for null data points to be plotted as white in an
%   IMAGESC plot, or similar.
%
%   [CBAR,CLIMS] = SMART_COLORBAR(CLIMS, CMAP_STYLE)
%
%   ARGUMENTS
%   * CLIMS: The desired colorbar limits
%   * CMAP_STYLE(Optional): Default jet(255). If a NX3 matrix compatible
%   with colormap, uses the colormap specified. If a string, analyzes the
%   string to find the colormap.
%
%   RETURNS
%   * CBAR: The new colobar with white background.
%   * CLIMS: The new color limits for the background. To make holes white,
%   set values to CLIMS(1)
%
% This code has been freely distributed by the authors. If used or
% modified, we would appreciate it if you cited our paper:
% Climer, J. R., Newman, E. L. and Hasselmo, M. E. (2013), Phase coding by 
%   grid cells in unconstrained environments: two-dimensional phase 
%   precession. European Journal of Neuroscience, 38: 2526–2541. doi: 
%   10.1111/ejn.12256
%
% RELEASE NOTES
%   v1.0 2014-10-15 Release (Jason Climer, jason.r.climer@gmail.com)
%   v0.1 2010-11-03 Originally written by Andrew Bogaard
%
% This file is part of pass_index.
%
% Copyright (c) 2014, Trustees of Boston University
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
%
% 1. Redistributions of source code must retain the above copyright notice, 
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright 
% notice, this list of conditions and the following disclaimer in the 
% documentation and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER 
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.   
if ~exist('cmap_style','var')
    cbar = jet(255);
else
    if ischar(cmap_style)
        eval(['cbar = ' cmap_style ';']); % initialize cbar
    else
        cbar = cmap_style;
    end
end

cbar = cat(1, [1 1 1], cbar);

nbins = size(cbar, 1)-1;

range = diff(clims)/(nbins-1)+diff(clims);

clims = [clims(2)-range, clims(2)];