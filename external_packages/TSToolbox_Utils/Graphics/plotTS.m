function axh = PlotTS(TSA, varargin)



opt_varargin = varargin;

defined_options = dictArray({ { 'MarkerSize', {5, {'numeric'}}}, 
                                              { 'Marker', {'k.', {'char'} } },
                                              { 'Offset', {0, {'numeric'} } },
                                              { 'Tlag', {0, {'numeric'} } },  
                                              { 'TUnits', {'s', {'char', 'TimeUnits'} } } 
                                              { 'YTicksOn', {1, {'numeric'} } }
                                              { 'FigureHandle', {[], {'numeric'} } } 
                                              { 'AxHandle', {[], {'numeric'} } } 
                                          } );
                                          
getOpt;

t = Range(TSA, TUnits)-Tlag;

if isempty(FigureHandle)
    FigureHandle = figure;
else
    figure(FigureHandle);
end

if ~isempty(AxHandle)
    axes(AxHandle);
end

plot(t, Offset * ones(size(t)), Marker, 'MarkerSize', MarkerSize);

if ~YTicksOn
    set(gca, 'YTick', []);
end

axh = gca;
