function [f, df] = ln_poisson_model(param,data,modelType)

%updated by LB to limit predictors to Pos,HD, and Theta

X = data{1}; % subset of A
Y = data{2}; % number of spikes

% compute the firing rate
u = X * param;
rate = exp(u);

% roughness regularizer weight - note: these are tuned using the sum of f,
% and thus have decreasing influence with increasing amounts of data
b_pos = 8e0; b_hd = 5e1; b_ego = 5e1; %b_theta = 5e1;

% % start computing the Hessian
% rX = bsxfun(@times,rate,X);
% hessian_glm = rX'*X;

%% find the P, H, S, or T parameters and compute their roughness penalties

% initialize parameter-relevant variables
J_pos = 0; J_pos_g = []; %J_pos_h = [];
J_hd = 0; J_hd_g = []; %J_hd_h = [];
J_ego = 0; J_ego_g = []; 
% J_theta = 0; J_theta_g = []; %J_spd_h = [];


% find the parameters
numPos = 400; numHD = 18; numEgo = 18; % hardcoded: number of parameters numSpd = 10; numTheta = 18; 
[param_pos,param_hd,param_ego] = find_param(param,modelType,numPos,numHD,numEgo); 

% compute the contribution for f, df, and the hessian
if ~isempty(param_pos)
    [J_pos,J_pos_g,~] = rough_penalty_2d(param_pos,b_pos); %J_pos_h
end

if ~isempty(param_hd)
    [J_hd,J_hd_g,~] = rough_penalty_1d_circ(param_hd,b_hd); % J_hd_h
end
% 
if ~isempty(param_ego)
    [J_ego,J_ego_g,~] = rough_penalty_1d_circ(param_ego,b_ego); % J_hd_h
end

%% compute f, the gradient, and the hessian

f = sum(rate-Y.*u) + J_pos + J_hd + J_ego;
df = real(X' * (rate - Y) + [J_pos_g; J_hd_g; J_ego_g]);
% hessian = hessian_glm + blkdiag(J_pos_h,J_hd_h,J_spd_h,J_theta_h);


%% smoothing functions called in the above script
function [J,J_g,J_h] = rough_penalty_2d(param,beta)

numParam = numel(param);
D1 = spdiags(ones(sqrt(numParam),1)*[-1 1],0:1,sqrt(numParam)-1,sqrt(numParam));
DD1 = D1'*D1;
M1 = kron(eye(sqrt(numParam)),DD1); M2 = kron(DD1,eye(sqrt(numParam)));
M = (M1 + M2);

J = beta*0.5*param'*M*param;
J_g = beta*M*param;
J_h = beta*M;

function [J,J_g,J_h] = rough_penalty_1d_circ(param,beta)

numParam = numel(param);
D1 = spdiags(ones(numParam,1)*[-1 1],0:1,numParam-1,numParam);
DD1 = D1'*D1;

% to correct the smoothing across first and last bin
DD1(1,:) = circshift(DD1(2,:),[0 -1]);
DD1(end,:) = circshift(DD1(end-1,:),[0 1]);

J = beta*0.5*param'*DD1*param;
J_g = beta*DD1*param;
J_h = beta*DD1;

function [J,J_g,J_h] = rough_penalty_1d(param,beta)

numParam = numel(param);
D1 = spdiags(ones(numParam,1)*[-1 1],0:1,numParam-1,numParam);
DD1 = D1'*D1;
J = beta*0.5*param'*DD1*param;
J_g = beta*DD1*param;
J_h = beta*DD1;

%% function to find the right parameters given the model type
function [param_pos,param_hd,param_ego] = find_param(param,modelType,numPos,numHD,numEgo)

param_pos = []; param_hd = []; param_ego = []; % param_dist = [];
   % One Factor
if all(modelType == [1 0 0])
    param_pos = param;
elseif all(modelType == [0 1 0])
    param_hd = param;
elseif all(modelType == [0 0 1])
    param_ego = param;
    %Two Factor
elseif all(modelType == [1 1 0 ])
    param_pos = param(1:numPos);
    param_hd = param(numPos+1:numPos+numHD);
elseif all(modelType == [1 0 1 ])
    param_pos = param(1:numPos);
    param_ego = param(numPos+1:numPos+numEgo);
elseif all(modelType == [0 1 1])
    param_hd = param(1:numHD);
    param_ego = param(numHD+1:numHD+numEgo);
    %Full Model
elseif all(modelType == [1 1 1]) 
    param_pos = param(1:numPos);
    param_hd = param(numPos+1:numPos+numHD);
    param_ego = param(numPos+numHD+1:numPos+numHD+numEgo);
end



