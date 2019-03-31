function cfaddbuttons(cffig)
%CFADDBUTTONS Add buttons to the curve fitting plot figure window

%   $Revision: 1.19.2.5 $  $Date: 2007/03/15 19:25:38 $
%   Copyright 2001-2007 The MathWorks, Inc.

% Clear out any old stuff
h0 = findall(cffig,'Type','uicontrol','Style','pushbutton');
if ~isempty(h0), delete(h0); end

% Define information for the buttons
strings = {xlate('Data...') xlate('Fitting...') xlate('Exclude...') xlate('Plotting...') xlate('Analysis...')};
tips = {xlate('Import, smooth, view, rename, and delete data') ...
        xlate('Add or change a fitted line') ...
        xlate('Create, view, and rename exclusion rules') ...
        xlate('Control which items are plotted') ...
        xlate('Predict, interpolate, differentiate, etc.')};
cbacks = {@cbkimport @cbkfit @cbkexclude @cbkplot @cbkanalyze};

% *** "tags" also defined in cfadjustlayout, and they must match
tags = {'cfimport' 'cffit' 'cfexclude' 'cfplot' 'cfanalyze'};

% Add the buttons to the figure
for j=1:length(strings)
   uicontrol(cffig,'Units','normalized', ...
                    'Position',[.2*j-.1,.9,.15,.05],...
                    'String',strings{j}, 'TooltipString',tips{j}, ...
                    'Callback',cbacks{j}, 'Tag',tags{j});
end

% ---------------------- callback for Import button
function cbkimport(varargin)
%CBKIMPORT Callback for Import button

delete(findall(gcbf,'Tag','cfstarthint'));

awtinvoke('com.mathworks.toolbox.curvefit.DataManager', 'showDataManager(I)', 0);

% ---------------------- callback for Exclude button
function cbkexclude(varargin)
%CBKEXCLUDE Callback for Exclude button

awtinvoke('com.mathworks.toolbox.curvefit.Exclusion', 'showExclusion');

% ---------------------- callback for Fit button
function cbkfit(varargin)
%CBKFIT Callback for Fit button

% createFitting creates a fitting panel only if one does not yet exist.
awtinvoke('com.mathworks.toolbox.curvefit.Fitting','showFitting');

% ---------------------- callback for Plot button
function cbkplot(varargin)
%CBKPLOT0 Callback for Plot button

awtinvoke('com.mathworks.toolbox.curvefit.Plotting', 'showPlotting');

% ---------------------- callback for Analysis button
function cbkanalyze(varargin)
%CBKANALYZE Callback for Analysis button

awtinvoke('com.mathworks.toolbox.curvefit.Analysis', 'showAnalysis');
