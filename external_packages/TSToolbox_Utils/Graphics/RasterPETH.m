function [fh, rasterAx, histAx] = RasterPETH(S, center, TStart, TEnd, varargin)

%  Display a PETH-raster plot for discrete
%  	
%  	USAGE:
%  	[fh, rasterAx, histAx] = RasterPETH(S, center, TStart, TEnd, 'ParameterName','Value')
%  	
%  	RasterPETH creates a raster plot, triggered on different event times.
%  	
%  	INPUTS:
%  	S - a ts (typically cell spikes)
%  	center - a ts of triggering times, represented as time 0 on the figure
%  	TStart - start time of the PETH (in sec), most often a negative value (to display spikes before center)
%  	TEnd   - end time of the PETH (in sec), most often a positive value (to display spikes before center)
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


font_name = 'Arial';
    font_size = 10;
    font_weight = 'bold';
    line_width = 2;

opt_varargin = varargin;

defined_options  = dictArray({ { 'RasterFraction', {0.7, {'numeric'}} }
                               { 'BinSize', {0.01, {'numeric'}}},
                                {'LineWidth', {1, {'numeric'} } },
                                {'Markers', { {}, {'cell'}} } ,
                                 {'MarkerTypes', { {}, {'cell'}}}, 
                                {'MarkerSize', { [], {'numeric'} } }
                                });
getOpt;

is = intervalSet(Range(center)+TStart, Range(center)+TEnd);

sweeps = intervalSplit(S, is, 'OffsetStart', TStart);
for iM = 1:length(Markers)
    Markers{iM} = (Range(Markers{iM}) - Range(center))/10; 
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
set(gca, 'XLim', [TStart TEnd]);
RasterPlot(sweeps, 'AxHandle', rasterAx, ...
    'FigureHandle', fh, ...
    'TStart', TStart, ...
    'TEnd', TEnd, ...
    'LineWidth', LineWidth, ...
    'Markers', Markers, ...
    'MarkerTypes', MarkerTypes, ...
    'MarkerSize', MarkerSize);

set(gca, 'Box', 'on');

axes(histAx);
set(histAx,'XTick',[]);

ss = oneSeries(sweeps);
sq = intervalRate(ss, regularIntervals(TStart, TEnd, BinSize));

dArea =  Data(sq)/length(sweeps);
%area(Range(sq, 's'), Data(sq)/length(sweeps), 'FaceColor', 'k');
bar(Range(sq, 's'), Data(sq)/length(sweeps),1, 'FaceColor', 'k');
set(gca, 'FontName', font_name);
set(gca, 'FontWeight', font_weight);
set(gca, 'FontSize', font_size);
set(gca, 'LineWidth', line_width);
set(gca, 'XLim', [TStart TEnd]);
if max(dArea) > 0
    set(gca, 'YLim', [0 max(dArea) * 1.2]);
end

yl = get(gca, 'YTick');
yl = yl(find(yl==floor(yl)));
set(gca, 'YTick', yl);

xl = get(rasterAx, 'XTick');
set(histAx,'XTick',xl);

fh = gcf;
