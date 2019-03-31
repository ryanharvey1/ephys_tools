function [fh, rasterAx, histAx] = RasterPETH(S, Pos, center, TStart, TEnd, varargin)



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
                                });
getOpt;

is = intervalSet(Range(center)+TStart, Range(center)+TEnd);

sweeps = intervalSplit(S, is, 'OffsetStart', TStart);
posSweeps = intervalSplit(pos, is, 'OffsetStart', TStart);
keyboard
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
set(gca, 'XLim', [TStart TEnd]/10);
RasterPosPlot(sweeps, 'AxHandle', rasterAx, ...
    'FigureHandle', fh, ...
    'TStart', TStart, ...
    'TEnd', TEnd, ...
    'LineWidth', LineWidth, ...
    'Markers', Markers, ...
    'MarkerTypes', MarkerTypes, ...
    'MarkerSize', MarkerSize);

set(gca, 'Box', 'on');
axes(histAx);
%  
%  ss = oneSeries(sweeps);
%  sq = intervalRate(ss, regular_interval(TStart, TEnd, BinSize));



dArea =  Data(sq)/length(sweeps);
area(Range(sq, 'ms'), Data(sq)/length(sweeps), 'FaceColor', 'k');
set(gca, 'FontName', font_name);
set(gca, 'FontWeight', font_weight);
set(gca, 'FontSize', font_size);
set(gca, 'LineWidth', line_width);
set(gca, 'XLim', [TStart TEnd]/10);
if max(dArea) > 0
    set(gca, 'YLim', [0 max(dArea) * 1.2]);
end
yl = get(gca, 'YTick');
yl = yl(find(yl==floor(yl)));
set(gca, 'YTick', yl);
fh = gcf;
