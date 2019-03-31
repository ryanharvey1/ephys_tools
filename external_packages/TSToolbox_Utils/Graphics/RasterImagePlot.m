function val = RasterPlot(S, varargin)

% RasterImagePlot2(S, height, bar, tstart, tend)
% 
% Plot a continuous function in a spike raster like form. 
% RasterImagePlot2 plots all waveforms in a matrix form.
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
  
% battaglia & peyrache 2001,2007
  

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
						{'GaussWidth', { 1, {'numeric'} } }

                                           });
                                       
getOpt;                                           

if isempty(FigureHandle)
    FigureHandle = figure;
else
    figure(FigureHandle);
end


if ~isempty(AxHandle)

    axes(AxHandle);
end

if isfinite(TStart )
  for i = 1:length(S)
    Si{i} = Restrict(S{i}, TStart, TEnd);
  end
else 
  Si = S;
end

l=length(Range(Si{1}));
lIx = 1;

for i = 2:length(Si)

	if l>length(Range(Si{i}));
		l = length(Range(Si{i}));
		lIx = i;
	end

end


matVal = zeros(length(Si),l);
  
for i = 1:length(Si)

  v = Data(Si{i});
 
  matVal(i,:) = v(1:l)';

end

gw = gausswin(GaussWidth);
gw = gw/sum(gw);
dt = median(diff(Range(Si{1})));

%  keyboard
matVal = convn(matVal',gw,'same')';
matVal = matVal(:,GaussWidth:end-GaussWidth);
times = Range(Si{lIx});
times = times(GaussWidth:end-GaussWidth);
val = tsd(times,matVal');

imagesc(times/10000,[1:length(Si)],matVal)
% xlim([TStart TEnd])
%  axis xy
hold on

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
        elseif length(mt)==1
            mr = mt;
            mc= 'r';
        else error('MarkerTypes elements should have one or two characters')
        end
        
        i = 1:length(Si);
        l = line(M, (i)*Height, 'LineStyle', 'none', 'Marker', mr, 'MarkerEdgeColor', mc);
        if ~isempty(MarkerSize)
            set(l, 'MarkerSize', 20);
        end
    end
end


