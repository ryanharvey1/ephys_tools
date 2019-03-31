function M = Rotate3dPlot(varargin)

% Rotate3dPlot  Rotates current 3D figure around the z axis
%
% M = Rotate3dPlot(varargin)
%
% INPUTS: 
%       varargin PARAMETERS:
%           rotStep (default 5 deg): how much to rotate by
%           maxRot (default 360 deg): how far to rotate
% OUTPUTS: 
%       M = movie of rotated plot.  
%
% ADR 1998, version v4.1, last modified 11/23/98 by ADR

% status PROMOTED
% v4.1 23 nov 98 fixed zooming axis problem


%-------------------
% parameters

rotStep = 5;
maxRot = 360;
Extract_varargin;

axisHandle = gca;
figHandle = ParentFigureHandle(gca);

nSteps = maxRot/rotStep;
set(gca, 'CameraViewAngleMode', 'manual');
if nargout == 1
   M = moviein(nSteps, figHandle);
   makingMovie = 1;
else
   makingMovie = 0;
end
   
%-------------------
% go
[az el] = view;
for iFrame = 1:nSteps
   az = az + rotStep;
   view(az, el);
   drawnow; 
   if makingMovie
      M(:,iFrame) = getframe(figHandle);
   end
end
