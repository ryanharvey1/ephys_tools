function [cftoolFit, results, isgood, hint] = cftoolcreatefit(...
    cftoolFit, panel, method, normalize, hint, dataset, outlier, fitname)
%CFTOOLCREATEFIT  Utility for CFTOOL to fit a curve to data
%
%   [CFTOOLFIT, RESULTS, ISGOOD, HINT] = CFTOOLCREATEFIT(...
%       CFTOOLFIT, PANEL, METHOD, NORMALIZE, HINT, DATASET, OUTLIER, FITNAME)

%   Copyright 2001-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.8 $    $Date: 2010/05/10 16:59:41 $ 

% Get the fit and dataset database objects
thefitdb = getfitdb;
thedsdb = getdsdb;

% convert the java UDDObject and set new properties.
cftoolFit=handle(cftoolFit);
cftoolFit.name=fitname;
cftoolFit.outlier=outlier;
cftoolFit.type=panel;

% Find the dataset in the database and set the dataset properties in the
% fit object
dataset = find( thedsdb, 'name', dataset );
cftoolFit.dataset = dataset.name;
cftoolFit.dshandle = dataset;

% pull the new fitoptions to use from where they have been stored on the fitdb
fitOptions=thefitdb.newOptions;

% set normalize
fitOptions.Normalize = normalize;

% set (or clear) the smoothing parameter
if isequal(method,'smoothing')
    sp = str2double( hint );
    if ~isnan(sp)
        fitOptions.smoothingparam=sp;
    else
        fitOptions.smoothingparam=[];
    end
end

% Set up other fitoptions
fitOptions.weights=dataset.weight;

% Exclude points if requested
if ~isequal(cftoolFit.outlier, cfGetNoneString() )
    outset = find(getoutlierdb,'name',cftoolFit.outlier);
    fitOptions.Exclude = cfcreateexcludevector(dataset,outset);
else
    fitOptions.Exclude = cfcreateexcludevector(dataset,[]);
end

% Need to ensure that the hint is setup properly so that it can be
% changed or interpreted.
cftoolFit.hint = hint;

% Need to ensure that the fit options in the cftool.fit object agree
% with the fit options we've just setup
cftoolFit.fitOptions = fitOptions;

%---------------------------------
% DO THE FIT
%---------------------------------
results = cftoolFit.doFit( method );

% Wrap the cftool.fit object up so the java code can handle it as a UDDObject.
isgood = cftoolFit.isGood;
hint = cftoolFit.hint;
cftoolFit=java(cftoolFit);

