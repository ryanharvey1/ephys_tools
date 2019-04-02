function  ClusterFunc_LimitSpikesByWaveforms(self)

% ClusterFunc_LimitSpikesByWaveforms (SpikelistCluster)
%
% Plots four subplots with waveforms, allows limits

%==========================================================
% PARAMETERS
%===========================================================
useSelfColor = false;

% GO

MCC = self.getAssociatedCutter();
MCC.StoreUndo('Limit by Waveforms');

MCS = MClust.GetSettings();
WV = self.GetWaveforms();
WVD = WV.data;
[nSpikes, nCh, nSamp] = size(WVD);

myLimits = {};
myUndo = MClustUtils.UndoSystem();

% if nCh~=4
%     error('MClust:Cutter', 'LimitSpikesByWaveforms assumes 4 channels.');
% end

F = figure('Name', ['Waveforms: ' self.name], ...
    'Tag', MCS.DeletableFigureTag, ...
    'Units','Normalized');

ax = tight_subplot(2,2,[.01 .03],[.1 .1],[.1 .1]); % added by ryan harvey

% ax = nan(nCh,1); 
% for iCh = 1:nCh
%     ax(iCh) = subplot(2,2,iCh);
%     set(ax(iCh), 'UserData', iCh);
% end
darkBackground(gcf,[0.2 0.2 0.2],[0.7 0.7 0.7])% Ryan Harvey addition (dark mode)

slider_nSpikesToDisplay = ...
    uicontrol('Parent', F, ...
    'Units','Normalized', ...
    'Position', [0 0 0.7 0.05], ...
    'Style', 'slider', ...
    'min', 1, 'max', nSpikes, ...
    'value', min(1000, nSpikes), ...
    'Callback', @(src,event)RedrawWaveforms);
DoneButton = ...
    uicontrol('Parent', F, ...
    'Units', 'Normalized', ...
    'Position', [0.85 0 0.15 0.05], ...
    'Style', 'PushButton', ...
    'String', 'DONE (save)', ...
    'Callback', @(src,event)Done);
CancelButton = ...
    uicontrol('Parent', F, ...
    'Units', 'Normalized', ...
    'Position', [0.70 0 0.15 0.05], ...
    'Style', 'PushButton', ...
    'String', 'CANCEL', ...
    'Callback', @(src,event)close(F));
LimitButton = uicontrol('Parent', F, ...
        'Units', 'Normalized', ...
        'Position', [0 0.95 0.15 0.05], ...
        'String', 'SingleLimit', ...
        'Callback', @(src,event)Limit);
UndoButton = ...
    uicontrol('Parent', F, ...
    'Units', 'Normalized', ...
    'Position', [0.85 0.95 0.15 0.05], ...
    'Style', 'PushButton', ...
    'String', 'Undo', ...
    'Callback', @(src,event)Undo);

% added by Ryan Harvey 2019
MultiLimitButton = uicontrol('Parent', F, ...
        'Units', 'Normalized', ...
        'Position', [.15 0.95 0.30 0.05], ...
        'String', 'MultiLimit (Enter when finished)', ...
        'Callback', @(src,event)MultiLimit);


RedrawWaveforms();

% added by Ryan Harvey 2019
% -------------------------------------
    function MultiLimit
        while true
            myUndo.StoreUndo(myLimits, 'limit');
            try
                [x,y] = ginput(2);
            catch
                return
            end
            if ~isempty(x) 
                jCh = get(gca, 'UserData');
                myLimits{end+1} = MClust.Limits.WaveformLimit(jCh, floor(x(1)), min(y), max(y));
                RedrawWaveforms();
                clear x y
            else
                return
            end
        end
    end
% -------------------------------------
    function Limit
        myUndo.StoreUndo(myLimits, 'limit');
        [x,y] = ginput(2);
        jCh = get(gca, 'UserData');
        myLimits{end+1} = MClust.Limits.WaveformLimit(jCh, floor(x(1)), min(y), max(y));
        RedrawWaveforms();
    end
% -------------------------------------
    function Done
        self.SetSpikes(ApplyLimits());
        close(F);
    end
% -------------------------------------
    function Undo
        myLimits = myUndo.PopUndo();
        RedrawWaveforms();
    end
% --------------------------------------
    function [S,iX] = ApplyLimits
        nLimits = length(myLimits);
        S = self.GetSpikes(); nS = length(S);
        keep = true(nLimits, nS);
        for iL = 1:nLimits
            [~,keep(iL,:)] = myLimits{iL}.ApplyLimit(WV, S);
        end
        keep = all(keep, 1);
        iX = 1:nS;
        S = S(keep); iX = iX(keep);
    end
% -------------------------------------
    function RedrawWaveforms

        [S, iX] = ApplyLimits();
        nS = length(S);
        nToPlot = min(nS, floor(get(slider_nSpikesToDisplay, 'Value')));
        r = randperm(nS, nToPlot);         
        for jCh = 1:nCh
            cla(ax(jCh)); hold on;
            if ~isempty(iX)
                h = plot(ax(jCh), 1:nSamp, squeeze(WVD(iX(r), jCh, :)),'w');
                if useSelfColor, set(h, 'color', self.color); end
            end
            darkBackground(gcf,[0.2 0.2 0.2],[0.7 0.7 0.7])% Ryan Harvey addition (dark mode)
            set(ax(jCh), 'YTick', [0], 'XTick', [], 'UserData', jCh);
        end
        for iL = 1:length(myLimits)
            plot(ax(myLimits{iL}.channel), ...
                myLimits{iL}.sample * [1 1], [myLimits{iL}.min myLimits{iL}.max], ...
                'color', 'k', 'LineWidth', 2);
        end
%         for jCh = 1:nCh
%             set(ax(jCh), 'YTick', [0], 'XTick', [], 'UserData', jCh);
% 		end
        axes(ax(1)); title(sprintf('showing %d of %d spikes', nToPlot, nS));
    end 
end