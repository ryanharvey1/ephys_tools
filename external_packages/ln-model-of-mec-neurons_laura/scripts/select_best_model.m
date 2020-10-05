%% Description
% This code will implement forward feature selection in order to determine
% the simplest model that best describes neural spiking. First, the
% highest-performing single-variable model is identified. Then, the
% highest-perfmoring double-variable model that includes the
% single-variable model is identified. This continues until the full model
% is identified. Next, statistical tests are applied to see if including
% extra variables significantly improves model performance. The first time
% that including variable does NOT signficantly improve performance, the
% procedure is stopped and the model at that point is recorded as the
% selected model.

% the model indexing scheme:
% ph, p, h

function selected_model = select_best_model(testFit)
testFit_mat = cell2mat(testFit);
LLH_values = reshape(testFit_mat(:,3),numFolds,numModels);

% find the best single model
singleModels = 2:3;
[~,top1] = max(nanmean(LLH_values(:,singleModels))); top1 = top1 + singleModels(1)-1;

% Single model is LLH1 & Full model is LLH2
LLH1 = LLH_values(:,top1); LLH2 = LLH_values(:,1);

[p_llh_12,~] = signrank(LLH2,LLH1,'tail','right');

if p_llh_12 < 0.05 % full model is sig. better
    selected_model = 1; % full model
else
    selected_model = top1; %single model
end

% re-set if selected model is not above baseline
pval_baseline = signrank(LLH_values(:,selected_model),[],'tail','right');

if pval_baseline > 0.05
    selected_model = NaN;
end
end

