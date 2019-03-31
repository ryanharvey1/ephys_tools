function [V,W] = FixTrackerJumps(X,Y, threshold_pixels, recursion_limit)
%
% [V,W] = FixTrackerJumps(X,Y, threshold_pixels, recursion_limit)
% 
% Remove regions (blocks) where the tracker obviously lost track of the headstage and 
% jumps > threshold_pixles pixels per frame.
% 
% This proceeds recursively up to  recursion_limit times. This way islands of up to recursion_limits frames can be removed.  
%  
%  inputs:
%  X,Y are(c)tsd objects returned by [X,Y]=LoadPosition('VT1.pascii');
% threshold_pixels .... distance in pixels between two consecutive frames above which the jump is regarded nonbiological
% recursion_limit  .... max number of times the algorithm is applied recursively (to remove islands of more than one frame)
%
%  outputs:
%  V,W are corrected tsd objects; the missing position coordinates removed, so frames are LOST!!
%
% PL 2001
% Version 0.2
% Note: Does not work following David's Video_Processing code! that code
% has its own algorithms for tracker jump removal
% ZN 2010 added figure ploting so could trouble shoot and visualize jumps


x = Data(X);
y = Data(Y);
ts = Range(X,'ts');

figure; plot(x,y, 'k-')
axis equal

next = 1;
while next > 0 
    ds = sqrt(diff(x).^2 + diff(y).^2);
    ix = find(ds > threshold_pixels) + 1;      % remove the second frame of a diff pair
    if isempty(ix)
        next = 0;
    else
        next = next + 1;
        disp([next, length(ix)]);
        hold on; plot(x(ix),y(ix),'r.')
        x(ix) = [];
        y(ix) = [];
        ts(ix) = [];
    end
            
    if next > recursion_limit
        next = 0;        
    end%if
end%while

V = tsd(ts,x);
W = tsd(ts,y);


