function r = cdot_hist(bins, freq, c, r, dotsize, ah, holdit)
% CDOT_HIST(bins, freq, c, r, dotsize, h, holdit)
% Creates a circular dot plot from an angular histogram
% 
% Input arguments:
% bins      - centre of histogram bins in degrees
% freq      - frequency of observations in each bin
% c         - dot color (default 'b')
% r         - radius of the histogram (by default at least 8*dotsize)
% dotsize   - diameter of each dot in pixels (default 10)
% ah        - handle to axes object (default gca)
% holdit    - optional argument. Set to '-a' to append this plot to the previous plotted one
%
% Output arguments:
% r         - radius of the final histogram (can be used as an input to circlabel) 
% 
% Example:
% bins      = 0:10:360;
% bearings  = randn(30, 1)*90+180;
% freq      = histcounts(bearings, bins);
% r         = cdot_hist(5:10:355, freq, 'k');
% cdot_label(0:30:330, r);
%
% To have the histogram automatically calculated, use circ_dotplot instead
% 
% See also cdot_label, cdot_plot

%% check input arguments and assign defaults
persistent oldfreq;

    if nargin<7 || isempty(holdit) || isempty(oldfreq)
        holdit = ''; 
        oldfreq = zeros(size(freq)); 
        fc = 'w'; % circle facecolor
    elseif strcmp(holdit, '-a')
        if length(bins)==length(oldfreq)
            fc = 'none';
        else
            error('Cannot append dot histogram if number of bins is different.')
        end
    else
        error('Unknown holdit option');
    end % there will be an error if oldfreq is not the same size as freq
    if nargin<6 || isempty(ah), ah = gca; end
    if nargin<5 || isempty(dotsize), dotsize = 10; end
    if nargin<4 || isempty(r)
        autosize = dotsize*round(1.5*max(freq)); 
        minsize  = dotsize*8;
        r        = max([autosize minsize]);
    end
    if nargin<3 || isempty(c), c = 'b'; end

%% Change angles so 0 is up top, and then transform into radians
bins = deg2rad(mod(450-bins, 360));

%% prepare axes
axes(ah); hold on; axis equal; axis off;

%% draw circle
rectangle('position', [-r -r r*2 r*2], 'curvature', [1 1], 'linewidth', 2, 'facecolor', fc);

%% draw individual dots
rds = 0.9*dotsize; % real dot size is 10% smaller to create a small gap between dots
for i = 1:length(bins) % for each bin
    offset = oldfreq(i);
    dc     = r - (0.5+offset:freq(i)+offset)*dotsize; % dot centre: centre position of each dot in this bin (as a distance from histogram centre)
    [x, y] = pol2cart(bins(i), dc); % calculate x/y coordinates from angle and radius
    
    for j = 1:length(x) % draw each dot
        rectangle('position', [x(j)-rds/2 y(j)-rds/2 rds rds], 'curvature', [1 1], 'edgecolor', c, 'facecolor', c);
    end
   
end

oldfreq = oldfreq + freq;