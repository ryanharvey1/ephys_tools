function n = getNumCells(ses)
% n = getNumCells(ses)
% 
% returns the number of cells for the given session object.

% The number of cells is obtained from the number of tfiles found:
% For nonsplit sessions:  n = number of tfiles
% For split sessions   :  n = number of tfiles/number of epochs
%
% PL Feb. 2003

if ses.split
    [n,ncols] = size(ses.tfilefullnames);
else
    n = length(ses.tfilefullnames)
end