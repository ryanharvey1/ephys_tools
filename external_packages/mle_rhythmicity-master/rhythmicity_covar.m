function [ params, confidence, stats, everything ] = rhythmicity_covar( data, session_duration, covars, varargin )
% RHYTHMICITY_COVAR Analyze rhythmicity as a linear function of covariates
%  Analyze rhythmicity allowing parameters (e.g. frequency) to vary
%  linearly with an arbitrarity number of covariates (e.g. speed).
% INPUT:
%   data: Can either be spike timestamps or a cell array of lags
%   session_duration: Duration of session
%   covars: Either a nXm matrix, where n is the number of spikes and m is
%       the number of covariates, or a cell array corresponding to each lag
%
% PARAMETERS:
%   In addition to the parameters in mle_rhythmicity:
%
%   plotit (true): If true, plots histogram and distribution estimates. Can
%       also be [false true], in which case we don't plot the static fit.
%   plot_axis (1): The axis (singular) to plot. To plot additional axes,
%       save output everything and see plot_rhythmicity_covar
%   $_covar (0): $ is replaced by a parameter name (a, tau, b, c, f, s). A
%       vector of the covariates that will modulate that parameter, with 0
%       corresponding with a constant and higher integers corresponding
%       with the columns of covars.
%   post_hocs (struct): struct with fields a, tau, b, c, f or s. Each field
%       contains the indices of the covariates for the leave-one-out
%       analysis. For leave multiple out analysis, use a cell array. For
%       example, to test f against columns 1 and 2+3, use
%       struct('f',{{1, [2 3]}}).
%   display ('verbose'): If 'verbose', prints information about the ongoing
%       analysis
%   covar_labels ({''}) - Labels for the covariate axis
%
% Copywrite (c) 2015,2016 Trustees of Boston University
% All rights reserved.
%
% This file is part of mle_rhythmicity revision 2.0. The last committed
% version of the previous revision is the SHA starting with 93862ac...
%
% This code has been freely distributed by the authors under the BSD
% license (http://opensource.org/licenses/BSD2-Clause). If used or
% modified, we would appreciate if you cited our paper:
%
% Climer JR, DiTullio R, Newman EL, Hasselmo ME, Eden UT. (2014),
% Examination of rhythmicity of extracellularly recorded neurons in the
% entorhinal cortex. Hippocampus, 25:460-473. doi: 10.1002/hipo.22383.

%% Parse input
PARAMS = {'tau','b','c','f','s','v','r'};% List of parameters - used for easier coding

display = [];

ip = inputParser;
warning('off','stats:mle:EvalLimit');

% mle_rhythmicity params
ip.addParamValue('max_lag', 0.6);
ip.addParamValue('epochs',[]);
ip.addParamValue('epochmode','lead');
ip.addParamValue('noskip',false);
ip.addParamValue('refractory_rise',false);
ip.addParamValue('f_range',[1 13]);
ip.addParamValue('alpha',0.05);
ip.addParamValue('t_axis',[ ]);
ip.addParamValue('runs',3);

% Covar specific params
ip.addParamValue('a_covar',NaN);
ip.addParamValue('post_hocs',struct);
ip.addParamValue('display','verbose');
ip.addParamValue('plotit',[true true]);
ip.addParamValue('plot_axis',1);
ip.addParamValue('covar_labels',{''});

% Pull out covariates
for i=1:numel(PARAMS)
    ip.addParamValue([PARAMS{i} '_covar'],0);
end

% Parse
ip.parse(varargin{:});

% Dump everything other than covariates to workspace
i = fields(ip.Results);
for i = i(cellfun(@(x)isempty(findstr('_covar',x)),i))'
    eval([i{1} '=ip.Results.' i{1} ';']);
end

% If noskip, drop s from PARAMS
if noskip
    PARAMS = PARAMS(~ismember(PARAMS,'s'));
end
if ~refractory_rise
    PARAMS = PARAMS(~ismember(PARAMS,'v'));
end

cov = struct;% Struct of covariate indicies

if isscalar(plotit)
    plotit = repmat(plotit,[1 2]);
end

% Load structure with covariate indicies
for i=1:numel(PARAMS)
    cov.(PARAMS{i}) = sort(unique(ip.Results.([PARAMS{i} '_covar'])));
    if ~isequal(cov.(PARAMS{i}),ip.Results.([PARAMS{i} '_covar']))
        warning([PARAMS{i} '_covar not unique/sorted. Output may not be shaped the same as the input.']);
    end
end
if ~isequal(ip.Results.a_covar,0)
    cov.r = ip.Results.a_covar;
end

for i=find(~ismember(PARAMS,fields(post_hocs)))
    post_hocs.(PARAMS{i}) = [];
end

for i=1:numel(PARAMS)
    if ~iscell(post_hocs.(PARAMS{i}))
        if isempty(post_hocs.(PARAMS{i})), post_hocs.(PARAMS{i}) = {};
        else post_hocs.(PARAMS{i}) = {post_hocs.(PARAMS{i})}; end;
    end
end

PARAMS = [PARAMS {'a'}];
for i=1:numel(PARAMS)
    for j=1:numel(post_hocs.(PARAMS{i}))
        if ~iscell(post_hocs.(PARAMS{i}))
           post_hocs.(PARAMS{i})= {post_hocs.(PARAMS{i})};
        end
        if ~isempty(post_hocs.(PARAMS{i}){j})&&~isequal(post_hocs.(PARAMS{i}){j},sort(unique(post_hocs.(PARAMS{i}){j})))
            post_hocs.(PARAMS{i}){j}=sort(unique(post_hocs.(PARAMS{i}){j}));
            warning(['post_hocs.' PARAMS{i} ' not unique/sorted. Output may not be the same shape as the input']);
        end        
    end
end
PARAMS = PARAMS(1:end-1);

if ismember('a',fields(post_hocs))
    post_hocs.r = post_hocs.a;
end


% Handle r-a conversion
if ~isnan(ip.Results.a_covar)
    r_covar = ip.Results.a_covar;
end

if ismember('a',fields(cov))
    cov.r = cov.a;
elseif isnan(cov.r)
    cov.r = 0;
end

if all(cellfun(@(x)isequal(cov.(x),0)||any(isnan(cov.(x))),PARAMS))
    exception = MException('mle_rhythmicity:rhythmicity_covar:NoCov',...
        ['No covariates defined in input. For non-covariate rhythmicity '...
        'analysis, please see mle_rhythmicity. See ']);
    throw(exception);
end

%%
if ~all(cellfun(@(x)ismember(0,cov.(x)),PARAMS))
    warning('Covarying parameters were used without a constant offset, if this was not intended it may results in bad fits and covariates cannot be shifted back from zero-mean.');
end

if isequal(display,'verbose')
    fprintf('\nrhythmicity_covar on %i events. Covariates are:',numel(data));
    for i=1:numel(PARAMS)
        fprintf('\n\t%s',PARAMS{i});
        if isequal(cov.(PARAMS{i}),0)
            fprintf(' is constant');
        else
            fprintf(' as a linear function of column');
            if sum(cov.(PARAMS{i})~=0)>1
                fprintf(['s ' strjoin(arrayfun(@num2str,cov.(PARAMS{i})(cov.(PARAMS{i})~=0),'UniformOutput',false),',')]);
            else
                fprintf(' %i',cov.(PARAMS{i})(cov.(PARAMS{i})~=0));
            end
            fprintf(' and ');
            if ~ismember(0,cov.(PARAMS{i}))
                fprintf('NO ');
            end
            fprintf('constant.');
        end
        
    end
    fprintf('\nRunning static fit...')
end

% Static fit
if any(plotit), HOLD = ishold; end

if plotit(1)% Static plotting
    if any(cellfun(@(x)isequal(x,'plotit'),varargin))
        varargin = varargin([1:find(cellfun(@(x)isequal(x,'plotit'),varargin))-1 find(cellfun(@(x)isequal(x,'plotit'),varargin))+2:end]);
    end
    pos = get(gca,'position');    
    if ~HOLD
        cla;
        axis off
    end
    if plotit(2)
        subplot('position',pos.*[1 1 19/44 1]);
    end
end

[~,~,stats0,everything0] = mle_rhythmicity(data,session_duration,'plotit',plotit,varargin{:});
if plotit(1)&&plotit(2)
    subplot('position',pos.*[1 1 19/44 1]+[pos(3)*25/44 0 0 0]);
end

% Pull out needed indices
inds = everything0.inds;
epochs = everything0.epochs;
lags = everything0.lags;
lags_list = everything0.lags_list;
n = numel(lags_list);

if isequal(display,'verbose')
    fprintf('done, %i lags, cell is ',numel(lags_list));
    if stats0.p_rhythmic>0.05
        fprintf('NOT ');
    end
    fprintf('rhythmic.');
end

%% Format covariates
% Align covars with lags
if ~iscell(covars)
    covars_lags =  cellfun(@(x)covars(x,:),inds,'UniformOutput',false);
else
    covars_lags = covars;
end

covars_list = cat(1,covars_lags{:});

% Shift the covariates to 0 mean - this makes the fit easier
shift = mean(covars_list);
covars_list = covars_list-repmat(shift,[size(covars_list,1) 1]);
covars_list = [ones(size(covars_list,1),1) covars_list];
shift = [0 shift];

% Log-likelihood wrapper function
LL_fun = @(cif_fun,cif_int,varargin)infbnd(sum(log(cif_fun(lags_list,varargin{:}))-log(cif_int(max_lag,varargin{:}))));

% Generate CIFs and CIF integral functions
if noskip
    if refractory_rise
        [cif_fun, cif_int] = cif_generator('noskip_rise');
    else
    [cif_fun, cif_int] = cif_generator('noskip');
    end
else
    if refractory_rise
         [cif_fun, cif_int] = cif_generator('full_rise');
    else
    [cif_fun, cif_int] = cif_generator('full');
    end
end

% Make fit bounds
lowerbound = [-inf 0 -inf 0 0];
upperbound = [inf 1 inf everything0.phat(5)*2 1];
if ~noskip
    lowerbound = [lowerbound 0];
    upperbound = [upperbound 1];
end
if refractory_rise
   lowerbound = [lowerbound(1:end-1) 0 lowerbound(end)];
   upperbound = [upperbound(1:end-1) inf upperbound(end)];
end
%%
% Make initial guess & large bounds
x0 = zeros(1,sum(cellfun(@(x)numel(cov.(x)),PARAMS)));

phat = everything0.phat;
phat = [phat(2:end) phat(1)/(1-phat(2))];
% clc
lbnd = -inf(1,numel(x0));
ubnd = inf(1,numel(x0));
for i=1:numel(PARAMS)
    lbnd(sum(arrayfun(@(j)numel(cov.(PARAMS{j})),1:i-1))+1)=lowerbound(i);
    ubnd(sum(arrayfun(@(j)numel(cov.(PARAMS{j})),1:i-1))+1)=upperbound(i);
end
%
for i=1:numel(PARAMS)
    x0(sum(arrayfun(@(j)numel(cov.(PARAMS{j})),1:i-1))+1)=phat(i);
end

%% Initial fits: allow covariates one parameter at a time
if isequal(display,'verbose')
    fprintf('\nFinding initial guess...');
end

for k=1:3% Iterate 3 times
    for i=find(~cellfun(@(x)isequal(cov.(x),0),PARAMS))% For each non-static parameter
        if isinf(upperbound(i))% The parameter is unbounded
            x0(sum(cellfun(@(x)numel(cov.(x)),PARAMS(1:i-1)))+(1:numel(cov.(PARAMS{i})))) = ... Set the initial value
                fminsearch(... Use fminsearch (unbounded)
                ... BEGIN FUN
                @(phat)...
                -passall(@(varargin)LL_fun(cif_fun,cif_int,varargin{:}),...log-likelihood function
                covar_wrapper(... Wrapper for parameters with covariates
                [x0(1:sum(cellfun(@(x)numel(cov.(x)),PARAMS(1:i-1))))... Keep other parameters fixed
                phat... The current covarying parameter
                x0(sum(cellfun(@(x)numel(cov.(x)),PARAMS(1:i)))+1:end)... Keep other parameters fixed
                ]...
                ,cov,covars_list)... END covar_wrapper
                )...
                ... END FUN
                , x0(sum(cellfun(@(x)numel(cov.(x)),PARAMS(1:i-1)))+(1:numel(cov.(PARAMS{i}))))... x0
                , optimset(... options
                'TolX',1e-10...
                ,'display','none'...
                ));
        else % Bounded covariates
            x0(sum(cellfun(@(x)numel(cov.(x)),PARAMS(1:i-1)))+(1:numel(cov.(PARAMS{i})))) = ... Set the initial value
                fmincon(... Use fmincon (bounded)
                ... BEGIN FUN
                @(phat)...
                -passall(@(varargin)LL_fun(cif_fun,cif_int,varargin{:}),...log-likelihood function
                covar_wrapper(...Wrapper for parameters with covariates
                [x0(1:sum(cellfun(@(x)numel(cov.(x)),PARAMS(1:i-1))))...Keep other parameters fixed
                phat...The current covarying parameter
                x0(sum(cellfun(@(x)numel(cov.(x)),PARAMS(1:i)))+1:end)...Keep other parameters fixed
                ]...
                ,cov,covars_list)...END covar_wrapper
                )...
                ... END FUN
                , x0(sum(cellfun(@(x)numel(cov.(x)),PARAMS(1:i-1)))+(1:numel(cov.(PARAMS{i}))))... x0
                , [covars_list(:,cov.(PARAMS{i})+1);-covars_list(:,cov.(PARAMS{i})+1)] ... A (linear constraints)
                , [ones(numel(lags_list),1)*upperbound(i);ones(numel(lags_list),1)*-lowerbound(i)] ... b (linear constraints)
                , [] ... Aeq
                , [] ... beq
                , [] ... ub
                , [] ... lb
                , [] ... nlcon
                , optimset(... options
                'algorithm','sqp'...
                ,'TolX',1e-10...
                ,'display','none'...
                ));
        end
    end
    
end

%% Final fits
% Make linear constraints for the final fits
A = ...
    arrayfun(@(k)... For each*
    [...
    zeros(2*numel(lags_list),sum(cellfun(@(x)numel(cov.(x)),PARAMS(1:k-1))))...% Zeros for irrelevant parameters
    [covars_list(:,cov.(PARAMS{i})+1);-covars_list(:,cov.(PARAMS{i})+1)]...% The covariates for this (k) bounded parameter (+ and -)
    zeros(2*numel(lags_list),sum(cellfun(@(x)numel(cov.(x)),PARAMS(k+1:end))))...% Zeros for irrelevant parameters
    ]...
    ,find(~cellfun(@(x)isequal(cov.(x),0),PARAMS)&~isinf(upperbound))...*constrained, covarying parameter
    ,'UniformOutput',false);
A = cat(1,A{:});
B = ...
    arrayfun(@(k)...For each*
    [ones(numel(lags_list),1)*upperbound(k);ones(numel(lags_list),1)*-lowerbound(k)]...% Bounds for this parameter (+upper, -lower)
    ,find(~cellfun(@(x)isequal(cov.(x),0),PARAMS)&~isinf(upperbound))...*constrained, covarying parameter
    ,'UniformOutput',false);
B = cat(1,B{:});

if isequal(display,'verbose')
    fprintf('done.\nFinal convergance...');
end

phat = fmincon(...
    ... BEGIN FUN
    @(phat)...
    -passall(@(varargin)LL_fun(cif_fun,cif_int,varargin{:})...log-likelihood function
    ,covar_wrapper(phat,cov,covars_list))...Wrapper for parameters with covariates
    ... END FUN
    ,x0 ... x0
    , A ... A
    , B ... B
    , [] ... Aeq
    , [] ... beq
    , lbnd ... lbnd
    , ubnd ... ubnd
    , [] ... nlcon
    , optimset(...
    'algorithm','interior-point'...
    ,'TolX',1e-10...
    ,'display','none'...'iter'...
    ));

LL = passall(@(varargin)LL_fun(cif_fun,cif_int,varargin{:})...
    ,covar_wrapper(phat,cov,covars_list));

if isequal(display,'verbose')
    fprintf('done.');
end

%% Post-hoc tests
post_hoc_results = struct;
if any(~cellfun(@(x)isempty(post_hocs.(x)),PARAMS))
    if isequal(display,'verbose')
        fprintf('\nRunning post hoc tests...');
    end
    for i=1:numel(PARAMS)
        for k=1:numel(post_hocs.(PARAMS{i}))
            j=post_hocs.(PARAMS{i}){k};
            if isequal(j,0)||~all(ismember(j,cov.(PARAMS{i})))
                warning('Covariates #%s not under %s, skipping.',strjoin(arrayfun(@num2str,j,'UniformOutput',false),', '),PARAMS{i});
            else
                post_hoc_results.(sprintf('%s%s',PARAMS{i},strjoin(arrayfun(@num2str,j,'UniformOutput',false),'_'))) = struct;
                cov2 = cov;
                cov2.(PARAMS{i}) = cov2.(PARAMS{i})(~ismember(cov2.(PARAMS{i}),j));
                post_hoc_results.(sprintf('%s%s',PARAMS{i},strjoin(arrayfun(@num2str,j,'UniformOutput',false),'_'))).phat = ...
                    mle(...
                    1,'logpdf',@(~,varargin)...
                    passall(@(varargin)LL_fun(cif_fun,cif_int,varargin{:})...
                    ,covar_wrapper([varargin{:}],cov2,covars_list))...
                    ,'start'...
                    ,phat([1:sum(cellfun(@(x)numel(cov.(x)),PARAMS(1:i-1)))...
                    sum(cellfun(@(x)numel(cov.(x)),PARAMS(1:i-1)))+find(~ismember(cov.(PARAMS{i}),j))...
                    (sum(cellfun(@(x)numel(cov.(x)),PARAMS(1:i)))+1):end]));
                post_hoc_results.(sprintf('%s%s',PARAMS{i},strjoin(arrayfun(@num2str,j,'UniformOutput',false),'_'))).LL = passall(@(varargin)LL_fun(cif_fun,cif_int,varargin{:})...
                    ,covar_wrapper(post_hoc_results.(sprintf('%s%s',PARAMS{i},strjoin(arrayfun(@num2str,j,'UniformOutput',false),'_'))).phat,cov2,covars_list));
                post_hoc_results.(sprintf('%s%s',PARAMS{i},strjoin(arrayfun(@num2str,j,'UniformOutput',false),'_'))).deviance = ...
                    2*(LL-...
                    post_hoc_results.(sprintf('%s%s',PARAMS{i},strjoin(arrayfun(@num2str,j,'UniformOutput',false),'_'))).LL);
                post_hoc_results.(sprintf('%s%s',PARAMS{i},strjoin(arrayfun(@num2str,j,'UniformOutput',false),'_'))).p = ...
                    1-chi2cdf(post_hoc_results.(sprintf('%s%s',PARAMS{i},strjoin(arrayfun(@num2str,j,'UniformOutput',false),'_'))).deviance,...
                    numel(j));
            end
        end
    end
end

%% Shift everything back
covars_list = covars_list+ones(size(lags_list))*shift;
if all(cellfun(@(x)ismember(0,cov.(x)),PARAMS(~cellfun(@(x)isequal(cov.(x),0),PARAMS))))
    for i=1:numel(PARAMS)
        phat(sum(cellfun(@(x)numel(cov.(x)),PARAMS(1:i-1)))+1) = ...
            phat(sum(cellfun(@(x)numel(cov.(x)),PARAMS(1:i-1)))+1) - ...
            sum(shift(cov.(PARAMS{i})+1).*...
            phat(sum(cellfun(@(x)numel(cov.(x)),PARAMS(1:i-1)))+(1:numel(cov.(PARAMS{i})))));
        for k=1:numel(post_hocs.(PARAMS{i}))
            j=post_hocs.(PARAMS{i}){k};
            cov2 = cov;
            cov2.(PARAMS{i}) = cov2.(PARAMS{i})(~ismember(cov.(PARAMS{i}),j));
            post_hoc_results.(sprintf('%s%s',PARAMS{i},strjoin(arrayfun(@num2str,j,'UniformOutput',false),'_'))).phat(...
                sum(cellfun(@(x)numel(cov2.(x)),PARAMS(1:i-1)))+1 ...
                ) = ...
                post_hoc_results.(sprintf('%s%s',PARAMS{i},strjoin(arrayfun(@num2str,j,'UniformOutput',false),'_'))).phat(...
                sum(cellfun(@(x)numel(cov2.(x)),PARAMS(1:i-1)))+1 ...
                ) - ...
                sum(shift(cov2.(PARAMS{i})+1).*...
                post_hoc_results.(sprintf('%s%s',PARAMS{i},strjoin(arrayfun(@num2str,j,'UniformOutput',false),'_'))).phat(...
                sum(cellfun(@(x)numel(cov2.(x)),PARAMS(1:i-1)))+(1:numel(cov2.(PARAMS{i})))));
        end
        
    end
end

if isequal(display,'verbose')
    fprintf('done.');
end

%% Packaging output
params = struct;
confidence = struct;
everything = struct;

for i=1:numel(PARAMS)% Pull out parameters going with each main parameter
    params.(PARAMS{i}) = ...
        phat(sum(cellfun(@(x)numel(cov.(x)),PARAMS(1:i-1)))+(1:numel(cov.(PARAMS{i}))));
end

if isequal(cov.b,0)% We can do a conversion
    PARAMS2 = [{'a'} PARAMS(~ismember(PARAMS,'r'))];
    params.a = (1-params.b)*params.r;% Add a to params
    
    % Do a conversions
    everything.phat = cellfun(@(x)params.(x),PARAMS2,'UniformOutput',false);
    everything.phat = cat(2,everything.phat{:});
    
    % Calculate a se
    everything.se = (diag(mlecov(phat,1,'logpdf',@(~,varargin)passall(@(varargin)...
        LL_fun(cif_fun,cif_int,varargin{:})...
        ,covar_wrapper([varargin{numel(cov.r)+1:end} [varargin{1:numel(cov.r)}]/(1-params.b)]...
        ,cov,covars_list)))))';
    
    PARAMS2 = [{'r'} PARAMS(~ismember(PARAMS,'r'))];
    for i=1:numel(PARAMS2)
        confidence.(PARAMS2{i}) = ...
            [-1;1]*everything.se(sum(cellfun(@(x)numel(cov.(x)),PARAMS2(1:i-1)))+(1:numel(cov.(PARAMS2{i}))))+repmat(everything.phat(sum(cellfun(@(x)numel(cov.(x)),PARAMS2(1:i-1)))+(1:numel(cov.(PARAMS2{i})))),[2 1]);
    end
else % b is static
    warning('Cannot do a conversion when b is not static.');
    everything.phat = phat;
    %    everything.se = diff(ci)/2/norminv(1-alpha);
    for i=1:numel(PARAMS)
        confidence.(PARAMS{i}) = ...
            [-1;1]*everything.se(sum(cellfun(@(x)numel(cov.(x)),PARAMS(1:i-1)))+(1:numel(cov.(PARAMS{i}))))+repmat(everything.phat(sum(cellfun(@(x)numel(cov.(x)),PARAMS(1:i-1)))+(1:numel(cov.(PARAMS{i})))),[2 1]);
    end
end

% Copy to everything
everything.ids = everything0.inds;
everything.epochs = everything0.epochs;
everything.phat0 = everything0.phat;
everything.lags = everything0.lags;
everything.lags_list = lags_list;
everything.covars_list = covars_list;
everything.max_lag = max_lag;
everything.noskip = noskip;
everything.refractory_rise = refractory_rise;
everything.session_duration = session_duration;
everything.LL0 = everything0.LL;
everything.cov = cov;
everything.shift = shift;
everything.shifted_back = ...
    all(cellfun(@(x)ismember(0,cov.(x)),PARAMS(~cellfun(@(x)isequal(cov.(x),0),PARAMS))));

% Calculate LL
everything.LL = LL;
everything.covar_labels = covar_labels;

% Do stats
stats.deviance = 2*(everything.LL-everything.LL0);
stats.df = numel(everything.phat)-numel(everything.phat0);
stats.p = 1-chi2cdf(stats.deviance,stats.df);
stats.post_hoc = post_hoc_results;
%% Plotit
if HOLD
    hold on;
else
    hold off;
end
if plotit(2)
    plot_rhythmicity_covar( everything, plot_axis, 'covar_label', covar_labels{plot_axis} );
end

end