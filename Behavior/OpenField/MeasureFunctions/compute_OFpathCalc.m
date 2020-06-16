function [ PathLength, IV, pathDist] = compute_OFpathCalc(x,y)
%compute_pathCalc calulates features of paths. 
% Input: 
%   x: x coordinates in cm 
%   y: y coordinates in cm 
%   'Frame_Rate' default 10 Hz
%   'Move_Threshold' default 3 cm/sec
%   'Smooth_Factor' default .8 - smoothing window size 80% of frame rate

%Output: 
%   PathLength: total distance in cm of path   
%   IV: Instantaneous velocity in cm/sec
%   pathDist: distance between consecutive points - will be size of xy
%   vectors - 1. 

% Note. For best results - input vectors should contain as few missing
% (NaN) values as possible. FixPos (ephys_tools\Utils\) function can be used to address missing
% values when needed. 


p = inputParser;
p.addParameter('Frame_Rate',10);
p.addParameter('Move_Threshold',3); %
p.parse(varargin{:});
fr = p.Results.Frame_Rate;
move = p.Results.Move_Threshold;

% First calulate distance between consecutive points
pathDist=sqrt(diff(x).^2 + diff(y).^2); % distance formula

% Get the velocity between consecutive points in cm per second
IV=pathDist*fr; % 10 = a frame rate of 10 hz      %instanteous velocity

% Compute the overall distance travelled 
PathLength=nansum(pathDist(pathIV>=move,1));%Total Path length for points greater than 3cm/s

end

