function setLegendable(aGraphic, isLegendable)
%SETLEGENDABLE   Allow or disallow a graphic from the legend.
%
%   SETLEGENDABLE(aGraphic, true) makes aGraphic appear on the legend for
%   its axes.
%
%   SETLEGENDABLE(aGraphic, false) removes aGraphic from the legend of its
%   axes.

%   Copyright 2011 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2011/04/11 16:09:28 $ 

anAnnotation = get( aGraphic, 'Annotation' );
aLegendEntry = get( anAnnotation', 'LegendInformation' );

if isLegendable
    set( aLegendEntry, 'IconDisplayStyle', 'on' );
else
    set( aLegendEntry, 'IconDisplayStyle', 'off' );
end

end