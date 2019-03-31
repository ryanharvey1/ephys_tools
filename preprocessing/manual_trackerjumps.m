function manual_trackerjumps(ts,x,y,StartofRec,EndofRec,path)
% manual_trackerjumps: Allows you to manually cut out xy coordinates that
% are outside the bounds of your maze. These can be caused by unplugs or
% if the rat jumps out. If you do not remove these points, your ratemap
% will be messed up. 
%
% Input:
%       ts
%       x
%       y
%       StartofRec: ts indicating the start points of your event
%       EndofRec: ts indicating the end points of your event
%
% This function depends on restrictMovement.m
%
% Ryan E Harvey (2018)

% 1 for elminating the points outside drawn shape & 0 for inside
restic_dir=1;

savets=[];
for i=1:length(StartofRec)
    % index out each event
    xtemp=x(ts>=StartofRec(i) & ts<=EndofRec(i));
    ytemp=y(ts>=StartofRec(i) & ts<=EndofRec(i));
    tstemp=ts(ts>=StartofRec(i) & ts<=EndofRec(i));
    % use the gui to cut out points
    [~,~,in]=restrictMovement(xtemp,ytemp,restic_dir);
    % save the ts where the tracker error exists
    savets=[savets,tstemp(in)];
end
% locate the index for each tracker error
in=ismember(ts,savets);
% save that index to your session folder so you won't have to do this again
% each time you run your data
save([path,filesep,'restrictxy.mat'],'in')
end

