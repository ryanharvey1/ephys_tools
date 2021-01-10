function [f, df] = ln_poisson_model_HDVel(param,data,modelType)

%updated by LB to limit predictors to Pos,HD, Linear Velocity, Angular
%Velocity

X = data{1}; % subset of A
Y = data{2}; % number of spikes

% compute the firing rate
u = X * param;
rate = exp(u);

% roughness regularizer weight - note: these are tuned using the sum of f,
% and thus have decreasing influence with increasing amounts of data
b_pos = 8e0; b_hd = 5e1; b_speed = 5e1; b_angvel = 5e1;

% % start computing the Hessian
% rX = bsxfun(@times,rate,X);
% hessian_glm = rX'*X;

%% find the P, H, S, or T parameters and compute their roughness penalties

% initialize parameter-relevant variables
J_pos = 0; J_pos_g = []; %J_pos_h = [];
J_hd = 0; J_hd_g = []; %J_hd_h = [];
J_speed = 0; J_speed_g = []; 
J_angvel = 0; J_angvel_g = []; %J_spd_h = [];


% find the parameters
numPos = 400; numHD = 18; numSpd = 10; numAngvel = 20; % hardcoded: number of parameters numSpd = 10; numTheta = 18; 
[param_pos,param_hd,param_speed,param_angvel] = find_param(param,modelType,numPos,numHD,numSpd, numAngvel); 

% compute the contribution for f, df, and the hessian
if ~isempty(param_pos)
    [J_pos,J_pos_g,~] = rough_penalty_2d(param_pos,b_pos); %J_pos_h
end

if ~isempty(param_hd)
    [J_hd,J_hd_g,~] = rough_penalty_1d_circ(param_hd,b_hd); % J_hd_h
end
% 
if ~isempty(param_speed)
    [J_speed,J_speed_g,~] = rough_penalty_1d(param_speed,b_speed); % J_hd_h
end
if ~isempty(param_angvel)
    [J_angvel,J_angvel_g,~] = rough_penalty_1d(param_angvel,b_angvel); % J_hd_h
end

%% compute f, the gradient, and the hessian

f = sum(rate-Y.*u) + J_pos + J_hd + J_speed + J_angvel;
df = real(X' * (rate - Y) + [J_pos_g; J_hd_g; J_speed_g; J_angvel_g]);
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
function [param_pos,param_hd,param_spd,param_angvel] = find_param(param,modelType,numPos,numHD,numSpd,numAngvel)

param_pos = []; param_hd = []; param_spd = []; param_angvel = [];
% Single parameter
if all(modelType == [1 0 0 0]) 
    param_pos = param;
elseif all(modelType == [0 1 0 0]) 
    param_hd = param;
elseif all(modelType == [0 0 1 0]) 
    param_spd = param;
elseif all(modelType == [0 0 0 1]) 
    param_angvel = param;
% Two Parameter
elseif all(modelType == [1 1 0 0])
    param_pos = param(1:numPos);
    param_hd = param(numPos+1:numPos+numHD);
elseif all(modelType == [1 0 1 0]) 
    param_pos = param(1:numPos);
    param_spd = param(numPos+1:numPos+numSpd);
elseif all(modelType == [1 0 0 1]) 
    param_pos = param(1:numPos);
    param_angvel = param(numPos+1:numPos+numAngvel);
elseif all(modelType == [0 1 1 0]) 
    param_hd = param(1:numHD);
    param_spd = param(numHD+1:numHD+numSpd);
elseif all(modelType == [0 1 0 1]) 
    param_hd = param(1:numHD);
    param_angvel = param(numHD+1:numHD+numAngvel);
elseif all(modelType == [0 0 1 1])  
    param_spd = param(1:numSpd);
    param_angvel = param(numSpd+1:numSpd+numAngvel);
% Three Parameter
elseif all(modelType == [1 1 1 0])
    param_pos = param(1:numPos);
    param_hd = param(numPos+1:numPos+numHD);
    param_spd = param(numPos+numHD+1:numPos+numHD+numSpd);
elseif all(modelType == [1 1 0 1]) 
    param_pos = param(1:numPos);
    param_hd = param(numPos+1:numPos+numHD);
    param_angvel = param(numPos+numHD+1:numPos+numHD+numAngvel);
elseif all(modelType == [1 0 1 1]) 
    param_pos = param(1:numPos);
    param_spd = param(numPos+1:numPos+numSpd);
    param_angvel = param(numPos+numSpd+1:numPos+numSpd+numAngvel);
elseif all(modelType == [0 1 1 1]) 
    param_hd = param(1:numHD);
    param_spd = param(numHD+1:numHD+numSpd);
    param_angvel = param(numHD+numSpd+1:numHD+numSpd+numAngvel);
% Full Model
elseif all(modelType == [1 1 1 1])
    param_pos = param(1:numPos);
    param_hd = param(numPos+1:numPos+numHD);
    param_spd = param(numPos+numHD+1:numPos+numHD+numSpd);
    param_angvel = param(numPos+numHD+numSpd+1:numPos+numHD+numSpd+numAngvel);
end