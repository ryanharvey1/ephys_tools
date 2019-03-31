function [fh, rasterAx, histAx, matVal] = ImagePETH(S, center, TStart, TEnd, varargin)

%  Displays a PETH/raster plot-like for continuous data.
%  	
%  	USAGE:
%  	[fh, rasterAx, histAx] = RasterPETH_c(S, center, TStart, TEnd, 'ParameterName','Value')
%  	
%  	RasterPETH_c creates a raster plot, triggered on different event times.
%  	
%  	INPUTS:
%  	S - a ts (typically cell spikes)
%  	center - a ts of triggering times, represented as time 0 on the figure
%  	TStart - start time of the PETH, most often a negative value (to display spikes before center)
%  	TEnd   - end time of the PETH, most often a positive value (to display spikes before center)
%  	
%  	OUTPUTS:
%  	fh - Figure handle
%  	rasterAx - Raster axis handle
%  	histAx - PETH axis handle
%  	
%  	OPTIONS:
%  	'RasterFraction' - fraction of the figure occupied by the raster (default 0.7)
%  	'BinSize' - bin size of the PETH, in timestamps (most often 10^-4 s., default 10)
%  	'LineWidth' - line width of the raster plot (default 1)
%  	'Markers' (optionnal) - array of ts corresponding to times at which markers are plotted.
%  	'MarkerTypes' (optionnal) - array of strings corresponding to Marker shape 
%  				    and color as used with plot function
%  	'MarkerSize' (optionnal) - vector of marker size
%  
%  	EXAMPLE:
%  	
%  	tsa = ts(10^6*rand(1000,1),'fixOrder',1)
%  	ct = [1:5:100]'*10^4;
%  	mk = ct+2000*rand(length(ct),1);
%  	ct = ts(ct);
%  	mk = ts(mk);
%  	
%  	figure(1),clf
%  	fh = RasterPETH(tsa, ct, -10000, 5000, 'Markers',{mk},'MarkerTypes',{'r*'})

% copyright (c) 2004 Francesco P. Battaglia, 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html


% Adrien Peyrache 2007, adapted from Francesco Battaglia codes


font_name = 'Arial';
    font_size = 10;
    font_weight = 'bold';
    line_width = 2;

opt_varargin = varargin;

defined_options  = dictArray({ { 'RasterFraction', {0.7, {'numeric'}} }
                               { 'BinSize', {10, {'numeric'}}},
                                {'LineWidth', {1, {'numeric'} } },
                                {'Markers', { {}, {'cell'}} } ,
                                 {'MarkerTypes', { {}, {'cell'}}}, 
                                {'MarkerSize', { [], {'numeric'} } }
                                {'GaussWidth', { 1, {'numeric'} } }
                                });
getOpt;

dt = median(diff(Range(S)));
st = TStart-GaussWidth*dt;
en = TEnd+GaussWidth*dt;

is = intervalSet(Range(center)+st, Range(center)+en);

sweeps = intervalSplit(S, is, 'OffsetStart', st);

for iM = 1:length(Markers)
    Markers{iM} = (Range(Markers{iM}) - Range(center))/10000; 
end

rf = RasterFraction * 0.8;
rasterAx = axes('position', [0.1 0.05 0.8 (rf+0.05)]);
histAx = axes('position', [0.1 (rf+0.15) 0.8 (0.75-rf)]);

fh = gcf;
axes(rasterAx);

set(gca, 'FontName', font_name);
set(gca, 'FontWeight', font_weight);
set(gca, 'FontSize', font_size);
set(gca, 'LineWidth', line_width);
matVal = RasterImagePlot(sweeps, 'AxHandle', rasterAx, ...
	'FigureHandle', fh, ...
	'TStart', st, ...
	'TEnd', en, ...
	'LineWidth', LineWidth, ...
	'Markers', Markers, ...
	'MarkerTypes', MarkerTypes, ...
	'MarkerSize', MarkerSize, ...
	'GaussWidth', GaussWidth);
set(gca, 'Box', 'on');
axes(histAx);


%  keyboard
dArea =  (mean(Data(matVal)'));
area(Range(matVal, 's'), dArea, 'FaceColor', 'k');
set(gca, 'FontName', font_name);
set(gca, 'FontWeight', font_weight);
set(gca, 'FontSize', font_size);
set(gca, 'LineWidth', line_width);
set(gca, 'XLim', [TStart TEnd]/10000);

yl = [min(0,1.2*min(dArea)) max(dArea) * 1.2];

if max(dArea) > 0
    set(gca, 'YLim', yl);
end

yM = max(yl);%floor(100*max(yl))/100;
ym = min(yl);%floor(100*min(yl))/100;
yl = [ym:(yM-ym)/3:yM];
set(gca, 'YTick', yl);
fh = gcf;
