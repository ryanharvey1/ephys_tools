function r = cdot_plot(bearings, binwidth, c, r, dotsize, ah, holdit)
% CDOT_PLOT(bearings, binwidth, c, r, dotsize, ah, holdit)
% Creates a circular dot plot from angular data
% 
% Input arguments:
% bearings  - angular data in degrees
% binwidth  - width of histogram bins (default: 10)
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
% bearings  = randn(30, 1)*90+180;
% r         = cdot_plot(bearings);
% cdot_label(0:30:330, r);
% 
% See also cdot_label, cdot_hist

%% check input arguments; if not defined, pass empty arguments to circ_dothist and defaults will be assigned there
if nargin<7, holdit = ''; end 
if nargin<6, ah = []; end
if nargin<5, dotsize = []; end
if nargin<4, r = []; end
if nargin<3, c = ''; end
if nargin<2 || isempty(binwidth), binwidth = 10; end
        
%% calculate histogram and plot 
bins      = 0:binwidth:360;
freq      = histcounts(bearings, bins);
r         = cdot_hist(bins(1:end-1)+binwidth/2, freq, c, r, dotsize, ah, holdit);