function [X,Y]=trackerjumps(x,y)
%
% [X,Y] = trackerjumps(X,Y)
% 
% Remove regions (blocks) where the tracker obviously lost track of the headstage and 
% jumps > threshold_pixles pixels per frame.
% 
% This proceeds recursively up to recursion_limit times. This way islands of up to recursion_limits frames can be removed.  
%  
%  inputs:
%  x,y position data
% 
%  outputs:
%  X,Y are corrected 
%
% Ryan Harvey 2018

recursion_limit=5;

X=x;
Y=y;
row=1:length(x);
rowidx=1:length(x);


next = 1;
while next > 0 
    ds = sqrt(diff(x).^2 + diff(y).^2);
    ix = find(ds > 15) + 1;      % remove the second frame of a diff pair
    if isempty(ix)
        next = 0;
    else
        next = next + 1;
        x(ix) = [];
        y(ix) = [];
        rowidx(ix)=[];
    end
            
    if next > recursion_limit
        next = 0;        
    end
end

% find index of removed frames
row(rowidx)=0;
row(row~=0)=1;
% 
% figure;plot(X,Y,'k');hold on
% plot(X(logical(row)),Y(logical(row)),'.r')

% make removed frames NaN for use in fillmissing below
X(logical(row))=NaN;
Y(logical(row))=NaN;

% fill in all tracker jumps
XY=fillmissing([X,Y],'pchip',1,'EndValues','nearest');
X=XY(:,1);
Y=XY(:,2);
% 
% figure;plot(X,Y,'k');hold on
% plot(X(logical(row)),Y(logical(row)),'.r')
end