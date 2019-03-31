function RasterPlot(PH, varargin)
% RasterPlot(S, height, bar, tstart, tend)
% 
% Plot spike trains as rasters.
% INPUT:
% S: tsdArray of ts objects, spike trains
% OPTIONS: 
% Height: the height allocated to each cell (includes spacing), default 1
% BarFraction: the fraction of height occupied by a cell, default 0.8
% TStart, TEnd: the beginning of end of time interval to be plot, default
% 'Markers': times denoting a special event for each trial which that will
% be marked by a symbol
% 'MarkerTypes': symobels to denote markers, (first character color
% (optional), second character symbol. Defaults to '*'
% the whole thing
  
% batta 2001 under construction
  

opt_varargin = varargin;

defined_options = dictArray({ {'Height', {1, {'numeric'} } }, 
                                              {'BarFraction', {0.8, {'numeric'} } },
                                               {'LineWidth', {1, {'numeric'} } },
                                              {'TStart', {NaN, {'numeric'} } },
                                              {'TEnd', {1e100, {'numeric'} } } 
                                              { 'FigureHandle', {[], {'numeric'} } } 
                                              { 'AxHandle', {[], {'numeric'} } } 
                                              {'Markers', {{}, {'cell'} } }
                                              {'MarkerTypes', { {}, {'cell' } } }
                                              {'MarkerSize', { [], {'numeric'} } }
					      {'MarkerEdgeColor', {{'r'}, {'cell'} } }
					      {'MarkerFaceColor', {{''}, {'cell'} } }
                                              {'TimeUnits', {{'ms'}, {'cell'} } }
                                              {'Offset', {0, {'numeric'} } }
                                           });

getOpt;
cmap = hsv;


if isempty(FigureHandle)
    FigureHandle = figure;
else
    figure(FigureHandle);
end


if ~isempty(AxHandle)

    axes(AxHandle);
end

if isfinite(TStart )
  for i = 1:length(PH)
    PHi{i} = Restrict(PH{i}, TStart, TEnd);
  end
else 
  PHi = PH;
end
for i = 1:length(PH)
  ph = mod(Data(PHi{i}),2*pi);
  COL{i} = floor((size(cmap,1)-1)*ph/(2*pi))+1;
  PHi{i} = Range(PHi{i});
end

for i = 1:length(PHi)
  rg = PHi{i};
  col = COL{i};
  for t=1:length(rg);
    x = [rg(t) rg(t)];
    y = [i i+Height]+Offset; 
    line(x, y, 'Color', cmap(col(t),:), 'LineWidth', LineWidth);
    hold on
  end
end



if ~isempty(Markers)
    if isempty(MarkerTypes)
        for iM  = 1:length(Markers)
            MarkerTypes{iM} = '*';
        end
    end
    
    for iM = 1:length(Markers)
        M = Markers{iM};
        mt = MarkerTypes{iM};
        if length(mt) == 2
            mr = mt(2);
            mc = mt(1);
	    mf = mt(1);
        elseif length(mt)==1
            mr = mt;
            mc= MarkerEdgeColor{1};
            mf= MarkerFaceColor{1};
        else error('MarkerTypes elements should have one or two characters')
        end
        
        i = 1:length(Si);
	if length(mf)
		l = line(M, (i+1/2)*Height, 'LineStyle', 'none', 'Marker', mr, 'MarkerEdgeColor', mc,'MarkerFaceColor',mf);
	else
		l = line(M, (i+1/2)*Height, 'LineStyle', 'none', 'Marker', mr, 'MarkerEdgeColor', mc);
	end
        if ~isempty(MarkerSize)
            set(l, 'MarkerSize', 20);
        end
    end
end


yl = get(gca, 'YTick');
yl = yl(find(yl==floor(yl)));
set(gca, 'YTick', yl);