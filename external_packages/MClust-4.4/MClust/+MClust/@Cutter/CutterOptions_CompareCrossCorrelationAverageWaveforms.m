function CutterOptions_CompareCrossCorrelationAverageWaveforms(self)
% plots cross correlations of average waveforms of only shows
%
% INPUTS
%
% OUTPUTS
%
% NONE

% Ryan Harvey 2019

%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%
OnePlot = true;

%%%%%%%%%%%%%%%%%%%%%%

MCS = MClust.GetSettings();
nClu = length(self.Clusters);
nToShow = sum(cellfun(@(C)~C.hide, self.Clusters));

iShow = 1;
[nR, nC] = MClustUtils.BestSubplots(nToShow);

stringsForLegend = {};

fig=figure('NumberTitle', 'off','Tag', MCS.DeletableFigureTag);
fig.Color=[.25 .25 .25];

for iC = 1:nClu
    if ~self.Clusters{iC}.hide
        
        WV = self.Clusters{iC}.GetWaveforms();
        [mWV, sWV, xr] = MClust.AverageWaveform(WV);
        
        waveforms(iC,:)=mWV(:);
    end
end


for i=1:size(waveforms,1)
    for ii=i:size(waveforms,1)
        wavecors(i,ii)=corr2(waveforms(i,:),waveforms(ii,:));
    end
end
wavecors(wavecors==0)=NaN;


h = heatmap(wavecors);
colormap jet
h.Title = 'WaveForm Correlations';
h.XLabel = 'clusters';
h.YLabel = 'clusters';
h.MissingDataLabel = [];
h.FontColor='w';
h.CellLabelFormat='%0.2g';
h.FontSize=10;

end