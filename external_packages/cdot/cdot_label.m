function cdot_label(angles, r, d, fontsize, ticklength, ah, labels)
% CDOT_LABEL(angles, r, d, fontsize, ticklength, ah, labels)
% Inserts ticks and labels into a circular plot
%
% Input arguments:
% angles     - which angles to mark (default [0 90 180 270])
% r          - radius of circular arena that should be labeled (default: radius of largest rectangle object in axes)
% d          - minimum distance between labels and arena (50% more for left and right sides) (default: 10% of r)
% fontsize   - label font size (default 10)
% ticklength - length of outwards ticks (default 2% of r, use negative values for inward ticks)
% ah         - handle axes object to plot into (default gca)
% labels     - cell array of labels; has to be the same size and order as angles
%
% Example:
% bearings  = randn(30, 1)*90+180;
% r         = cdot_plot(bearings);
% cdot_label(0:30:330, r);
%
% See also cdot_hist, cdot_plot

%% Calculate appropriate fontsize and distances from circle
if nargin<7, labels={}; end
if nargin<6 || isempty(ah),         ah = gca; end
if nargin<2 || isempty(r),          r = nested_autodetect_size; end
if nargin<5 || isempty(ticklength), ticklength = r*0.02; end
if nargin<4 || isempty(fontsize),   fontsize = 10; end
if nargin<3 || isempty(d),          d = 0.1*r; end
if nargin<1 || isempty(angles),     angles = 0:90:270; end

for i=1:length(angles)
    if isempty(labels)
        nested_label(angles(i), [num2str(angles(i)) '\circ']);
    else
        nested_label(angles(i), labels{i});
    end
    nested_tick(angles(i));
end

%% nested functions
    function nested_label(theta, lab)
        % theta is 0 for up, 90 for right
        t = deg2rad(mod(450-theta, 360));
        [x, y] = pol2cart(t, r+(1+0.5*cos(t)^2)*d);
        text(x, y, lab, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle', 'FontUnits', 'points', 'FontSize', fontsize, 'parent', ah);
    end

    function nested_tick(theta)
        t = deg2rad(mod(450-theta, 360));
        polar([t t], [r r+ticklength], 'k');
    end

    function auto_r = nested_autodetect_size
        allrects = findobj('parent', ah, 'type', 'rectangle');
        if isempty(allrects)
            error('Size auto-detection failed, no rectangle objects found in these axes');
        end
        p = get(allrects, 'Position');
        diameters = cellfun(@(x) max(x(3:4)), p);
        auto_r = max(diameters)/2;
    end
end %main function