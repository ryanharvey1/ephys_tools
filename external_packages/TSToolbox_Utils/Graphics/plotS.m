function axh = PlotTS(S, varargin)

opt_varargin = varargin;

defined_options = dictArray({ { 'MarkerSize', {5, {'numeric'}}}, 
                                              { 'Marker', {'k.', {'char'} } },
                                              { 'OffsetShift', {1, {'numeric'} } },
                                              { 'Offset', {0, {'numeric'} } },
                                              { 'Tlag', {0, {'numeric'} } },  
                                              { 'TUnits', {'s', {'char', 'TimeUnits'} } } 
                                              { 'YTicksOn', {1, {'numeric'} } }
                                              { 'FigureHandle', {[], {'numeric'} } } 
                                              { 'AxHandle', {[], {'numeric'} } } 
                                          } );
getOpt;


if isempty(FigureHandle)
    FigureHandle = figure;
    hold on
else
    figure(FigureHandle)
end

for i=1:length(S)
	
	plotTS(S{i},'MarkerSize',MarkerSize,'Marker',Marker,'Offset',OffsetShift*i+Offset,...
		'Tlag',Tlag,'TUnits',TUnits,'YTicksOn',YTicksOn,...
		'FigureHandle',FigureHandle,'AxHandle',AxHandle);
	hold on
end