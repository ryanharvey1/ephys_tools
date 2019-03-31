function [C,B,CC,db] = MCrossCorr(t1,t2,dx,nbins)
% function to check Francescos C CrossCorr
%
%  t1,t2 ... input timestamp vectors of cell 1 and 2 (in 0.1 msec units)
%  dx    ... lag binsize in msec
%  nbins ... number of lag bins (should be odd - if nbins is even, nbins+1 C-values are calculated and returned) 
%            
%
%  C ...    cross-correlation or PETH of t2 around events given at t1 (t1 = events, t2 = snapshots) 
%  B ...    array of timestamps of lower-edges of - 
%  CC ...   unnormalized Coincidence Counts
%  db.DT2 ...   t2(i2)-TL, the T difference of the first spike i2 of t2 inside the [TL,TR] window (for debugging) 
%  db.DT1 ...   t2(i2-1)-TL, the T difference of the the spike PRECEDING the first spike i2 of t2 inside the [TL,TR] window (for debugging) 
%

% nbins needs to be odd
if 2*floor(nbins/2)== nbins
    nbins = nbins+1; 
end
nt1 = length(t1);   % # of spikes in t1 vector
nt2 = length(t2);   % # of spikes in t2 vector
W = dx*nbins/2;     % half-width of gliding obersvation window
CC = zeros(nbins,1); % initialize Coincidence Counter C
db.DT1 = zeros(nt1,1);  % for debugging
db.DT2 = zeros(nt1,1);  % for debugging

if nargout > 1
    B = ((-W+dx/2):dx:(W-dx/2))';   % array of bin-centers in 0.1 msec units
end

dx = dx*10;      % goto 0.1 msec units (same as t1 and t2)
W = W*10;        % recompute the window size in 0.1 msec units

i2 = 1;
for i1=1:nt1     % loop over spikes i1 of t1 vector
    TL = t1(i1) - W;    % shift observation window 2W = [TL,TR] so that it is centered at spike i1; 
                        % the obersvation window has length 2W and is now the open interval 2W = [t1(i1)-W, t1(i1)+W[  
    while t2(i2) < TL && i2 < nt2
        i2 = i2 + 1;    % if spike i2 is left of TL, search forward for the first spike i2 in t2 that is within the 2W window
    end
    while i2 > 1 && t2(i2-1) > TL 
        i2 = i2 - 1;    % if the PREVIOUS spike i2-1 is right of TL (inside 2W), search backwards for 
                        % first PRECEDING spike i2-1 in t2 that is to the left of TL. The following spike should 
                        % now be the first one within 2W 
    end
    
    % compute debugging info
    db.DT2(i1) = t2(i2) - TL;
    if i2 > 1
        db.DT1(i1) = t2(i2-1) - TL;
    else
        db.DT1(i1) = NaN;
    end
        
    % at this point i2 should be the first spike inside the 2W window no
    % matter if it came from a forward or backward search
    
    % fill all spikes i2,i2+1, ....  in a histogram of nbins bins (with binsize dx) starting
    % at TL and ending at TR = TL + nbins*dx
    TR = TL;        % initial condition where the moving TR equals TL
    l = i2;         % initialize temporary i2 index l for filling the histogram
    for j=1:nbins   % loop over bins j
        k = 0;      % initialize spike counter k in bin j
        TR = TR + dx;   % set the upper edge of the current bin j
        while l < nt2 && t2(l) < TR   % count the next l spikes as long as they are < TR
            k = k+1;    % increment counter k
            l = l+1;    % goto next spike l+1
        end
        CC(j) = CC(j) + k;
    end % j
    
end %i1

% normalization
C = 10000*CC/(dx*nt1);     % C is normalized such that it represents a PETH of spikes t2 centered at events t1. 
                           % The ordinate is then the avg. firing rate in Hz


%   /* cross correlations */
%   
%   w = ( (double)nbins / 2) * binsize;
% 
%   for(i1 = 0; i1 < nt1; i1++)
%   {
%       lbound = t1[i1] - w;
%       while(t2[i2] < lbound && i2 < nt2 -1)
%         i2++;
%       while(t2[i2-1] > lbound && i2 > 1)
%         i2--;
%       rbound = lbound;
%       l = i2;
%       
%       for(j = 0; j < nbins; j++)
%       {
%           k = 0;
%           rbound += binsize;
%           while(t2[l] < rbound && l < nt2-1)
%           {
%               l++;
%               k++;
%           }
%           C[j] += k;
%       }
%   }
%     
%   for(j = 0; j < nbins; j++)
%     C[j] /= nt1 * binsize / 10000;
