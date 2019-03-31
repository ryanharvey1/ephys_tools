function [peak, base, phase, gof] = cosfit8(x,orig,peaktopeak,peakoff,zeroline,dirt)

%COSFIT8    Create plot of datasets and fits
%   COSFIT8(X,ORIG,DIRT)
%   Creates a plot, similar to the plot in the main curve fitting
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with cftool
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of datasets:  1
%   Number of fits:  1

 
% Data from dataset "orig vs. x with dirt":
%    X = x:
%    Y = orig:
%    Weights = dirt:
%
% This function was automatically generated on 10-Jun-2009 19:25:01

% Set up figure to receive datasets and fits
f_ = clf;
figure(f_);
set(f_,'Units','Pixels','Position',[1772 364 680 484]);
legh_ = []; legt_ = {};   % handles and text for legend
xlim_ = [Inf -Inf];       % limits of x axis
ax_ = axes;
set(ax_,'Units','normalized','OuterPosition',[0 0 1 1]);
set(ax_,'Box','on');
axes(ax_); hold on;

 
% --- Plot data originally in dataset "orig vs. x with dirt"
x = x(:);
orig = orig(:);
dirt = dirt(:);
h_ = line(x,orig,'Parent',ax_,'Color',[0.333333 0 0.666667],...
     'LineStyle','-', 'LineWidth',1,...
     'Marker','.', 'MarkerSize',12);
xlim_(1) = min(xlim_(1),min(x));
xlim_(2) = max(xlim_(2),max(x));
legh_(end+1) = h_;
legt_{end+1} = 'orig vs. x with dirt';

% Nudge axis limits beyond data limits
if all(isfinite(xlim_))
   xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
   set(ax_,'XLim',xlim_)
end


% --- Create fit "fit 1"
ok_ = ~(isnan(x) | isnan(orig) | isnan(dirt));
st_ = [peaktopeak peakoff zeroline ];
ft_ = fittype('a*cos(x+b)+c' ,...
     'dependent',{'y'},'independent',{'x'},...
     'coefficients',{'a', 'b', 'c'});

% Fit this model using new data
[cf_, gof] = fit(x(ok_),orig(ok_),ft_ ,'Startpoint',st_,'Weight',dirt(ok_),'Lower',[0 -Inf 0], 'Upper', [Inf Inf 20]);
peak=cf_.a;
phase=cf_.b;
base=cf_.c;

% Or use coefficients from the original fit:
if 0
   cv_ = {-0.2906144055756, 2.924742804905, 7.52847430406};
   cf_ = cfit(ft_,cv_{:});
end

% Plot this fit
h_ = plot(cf_,'fit',0.95);
legend off;  % turn off legend from plot method call
set(h_(1),'Color',[1 0 0],...
     'LineStyle','-', 'LineWidth',2,...
     'Marker','none', 'MarkerSize',6);
legh_(end+1) = h_(1);
legt_{end+1} = 'fit 1';

% Done plotting data and fits.  Now finish up loose ends.
hold off;
h_ = legend(ax_,legh_,legt_,'Location','NorthEast');  
set(h_,'Interpreter','none');
ylabel(ax_,'');               % remove x label
%lambda=(.5/abs(peak))*nanmean(dirt(1:8,2));
xlabel(ax_,['r2=' num2str(round(gof.rsquare*1000)/1000) '; pd=' num2str(round((2*pi-phase)*100)/100) '; bl=' num2str(round(base*100)/100)]); % remove y label
% title(['preferred direction = ' num2str(phi)])






