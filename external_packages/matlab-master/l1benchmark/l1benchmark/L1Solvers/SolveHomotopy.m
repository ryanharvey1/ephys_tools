%% This function is modified from Matlab Package: L1-Homotopy

% BPDN_homotopy_function.m
% 
% Solves the following basis pursuit denoising (BPDN) problem
% min_x  \lambda ||x||_1 + 1/2*||b-Ax||_2^2
%
% Inputs: 
% A - m x n measurement matrix
% b - measurement vector
% lambda - final value of regularization parameter
% maxiter - maximum number of homotopy iterations
%
% Outputs:
% x_out - output for BPDN
% total_iter - number of homotopy iterations taken by the solver
% 
% Written by: Salman Asif, Georgia Tech
% Email: sasif@ece.gatech.edu
%
%-------------------------------------------+
% Copyright (c) 2007.  Muhammad Salman Asif 
%-------------------------------------------+

function [x_out, total_iter, timeSteps, errorSteps, epsSteps] = SolveHomotopy(A, b, varargin)

global N gamma_x z_x  xk_temp del_x_vec pk_temp dk epsilon isNonnegative

t0 = tic ;

lambda = 1e-6;
maxiter = 100;
isNonnegative = false;
verbose = false;
xk_1 = [];

STOPPING_TIME = -2;
STOPPING_GROUND_TRUTH = -1;
STOPPING_DUALITY_GAP = 1;
STOPPING_SPARSE_SUPPORT = 2;
STOPPING_OBJECTIVE_VALUE = 3;
STOPPING_SUBGRADIENT = 4;
STOPPING_DEFAULT = STOPPING_OBJECTIVE_VALUE;

stoppingCriterion = STOPPING_DEFAULT;

% Parse the optional inputs.
if (mod(length(varargin), 2) ~= 0 ),
    error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
end
parameterCount = length(varargin)/2;

for parameterIndex = 1:parameterCount,
    parameterName = varargin{parameterIndex*2 - 1};
    parameterValue = varargin{parameterIndex*2};
    switch lower(parameterName)
        case 'stoppingcriterion'
            stoppingCriterion = parameterValue;
        case 'initialization'
            xk_1 = parameterValue;
            if ~all(size(xk_1)==[n,1])
                error('The dimension of the initial x0 does not match.');
            end
        case 'groundtruth'
            xG = parameterValue;
        case 'lambda'
            lambda = parameterValue;
        case 'maxiteration'
            maxiter = parameterValue;
        case 'isnonnegative'
            isNonnegative = parameterValue;
        case 'tolerance'
            tolerance = parameterValue;
        case 'verbose'
            verbose = parameterValue;
        case 'maxtime'
            maxTime = parameterValue;
        otherwise
            error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
    end
end
clear varargin

timeSteps = nan(1,maxiter) ;
errorSteps = nan(1,maxiter) ;
epsSteps = nan(1,maxiter);

[K, N] = size(A);

% Initialization of primal and dual sign and support
z_x = zeros(N,1);
gamma_x = [];       % Primal support

% Initial step
Primal_constrk = -A'*b;

if isNonnegative
    [c i]  = min(Primal_constrk);
    c = max(-c, 0);
else
    [c i] = max(abs(Primal_constrk));
end

epsilon = c;
nz_x = zeros(N,1);
if isempty(xk_1)
    xk_1 = zeros(N,1);
    gamma_xk = i;
else
    gamma_xk = find(abs(xk_1)>eps*10);
    nz_x(gamma_xk) = 1;
end
f = epsilon*norm(xk_1,1) + 1/2*norm(b-A*xk_1)^2;
z_x(gamma_xk) = -sign(Primal_constrk(gamma_xk));
%Primal_constrk(gamma_xk) = sign(Primal_constrk(gamma_xk))*epsilon;

z_xk = z_x;

% loop parameters
iter = 0;
out_x = [];
old_delta = 0;
count_delta_stop = 0;

AtgxAgx = A(:,gamma_xk)'*A(:,gamma_xk);
iAtgxAgx = inv(A(:,gamma_xk)'*A(:,gamma_xk));
while iter < maxiter
    iter = iter+1;
    % warning('off','MATLAB:divideByZero')

    gamma_x = gamma_xk;
    z_x = z_xk;
    x_k = xk_1;

    %%%%%%%%%%%%%%%%%%%%%
    %%%% update on x %%%%
    %%%%%%%%%%%%%%%%%%%%%
    
    % Update direction
    del_x = iAtgxAgx*z_x(gamma_x);
    del_x_vec = zeros(N,1);
    del_x_vec(gamma_x) = del_x;

    %dk = A'*(A*del_x_vec);
    Asupported = A(:,gamma_x);
    Agdelx = Asupported*del_x;
    dk = A'*Agdelx;
    
    %%% CONTROL THE MACHINE PRECISION ERROR AT EVERY OPERATION: LIKE BELOW. 
    pk_temp = Primal_constrk;
    gammaL_temp = find(abs(abs(Primal_constrk)-epsilon)<min(epsilon,2*eps));
    pk_temp(gammaL_temp) = sign(Primal_constrk(gammaL_temp))*epsilon;
    
    xk_temp = x_k;
    xk_temp(abs(x_k)<2*eps) = 0;
    %%%---
    
    % Compute the step size
    [i_delta, delta, out_x] = update_primal(out_x);
    
    if old_delta < 4*eps && delta < 4*eps
        count_delta_stop = count_delta_stop + 1;
        
        if count_delta_stop >= 500
            if verbose
                disp('stuck in some corner');
            end
            break;
        end
    else
        count_delta_stop = 0;
    end
    old_delta = delta;
    
    xk_1 = x_k+delta*del_x_vec;
    Primal_constrk = Primal_constrk+delta*dk;
    epsilon_old = epsilon;
    epsilon = epsilon-delta;

    if epsilon <= lambda;
%         xk_1 = x_k + (epsilon_old-lambda)*del_x_vec;
%         disp('Reach prescribed lambda in SolveHomotopy.');

        gamma_x0 = find(abs(xk_1) > 1e-9);
        AtgxAgx0 = A(:,gamma_x0)'*A(:,gamma_x0);
        x_temp = AtgxAgx0 \ (A(:,gamma_x0)' * b);
        xk_1 = zeros(N,1);
        xk_1(gamma_x0) = x_temp;

        timeSteps(iter) = toc(t0) ;
        errorSteps(iter) = norm(xk_1-xG) ;
        epsSteps(iter) = epsilon;       
        break;
    end
    
    timeSteps(iter) = toc(t0) ;
    errorSteps(iter) = norm(xk_1-xG) ;
    epsSteps(iter) = epsilon;
    
    % compute stopping criteria and test for termination
    keep_going = true;
    
    switch stoppingCriterion
        case STOPPING_GROUND_TRUTH
            keep_going = norm(xk_1-xG)>tolerance;
        case STOPPING_SPARSE_SUPPORT
            if delta~=0
                nz_x_prev = nz_x;
                nz_x = (abs(xk_1)>eps*10);
                num_nz_x = sum(nz_x(:));
                num_changes_active = (sum(nz_x(:)~=nz_x_prev(:)));
                if num_nz_x >= 1
                    criterionActiveSet = num_changes_active / num_nz_x;
                    keep_going = (criterionActiveSet > tolerance);
                end
            end
        case STOPPING_DUALITY_GAP
            error('Duality gap is not a valid stopping criterion for Homotopy.');
        case STOPPING_OBJECTIVE_VALUE
            if delta~=0
                % continue if not yeat reached target value tolA
                prev_f = f;
                f = lambda*norm(xk_1,1) + 1/2*norm(b-Asupported*xk_1(gamma_x))^2;
                keep_going = (abs((prev_f-f)/prev_f) > tolerance);
            end
        case STOPPING_SUBGRADIENT
            keep_going = norm(delta*del_x_vec)>tolerance;
        case STOPPING_TIME
            keep_going = timeSteps(iter) < maxTime ;
        otherwise,
            error('Undefined stopping criterion');
    end % end of the stopping criteria switch

    
%     if keep_going && norm(xk_1 - x_k)<100*eps
%         if verbose
%             disp('The iteration is stuck.');
%         end
%         keep_going = false;
%     end
    
    if ~keep_going
        break;
    end

    if ~isempty(out_x)
        % If an element is removed from gamma_x
        len_gamma = length(gamma_x);
        outx_index = find(gamma_x==out_x(1));
        gamma_x(outx_index) = gamma_x(len_gamma);
        gamma_x(len_gamma) = out_x(1);
        gamma_x = gamma_x(1:len_gamma-1);
        gamma_xk = gamma_x;
        
        rowi = outx_index; % ith row of A is swapped with last row (out_x)
        colj = outx_index; % jth column of A is swapped with last column (out_lambda)
        AtgxAgx_ij = AtgxAgx;
        temp_row = AtgxAgx_ij(rowi,:);
        AtgxAgx_ij(rowi,:) = AtgxAgx_ij(len_gamma,:);
        AtgxAgx_ij(len_gamma,:) = temp_row;
        temp_col = AtgxAgx_ij(:,colj);
        AtgxAgx_ij(:,colj) = AtgxAgx_ij(:,len_gamma);
        AtgxAgx_ij(:,len_gamma) = temp_col;
        iAtgxAgx_ij = iAtgxAgx;
        temp_row = iAtgxAgx_ij(colj,:);
        iAtgxAgx_ij(colj,:) = iAtgxAgx_ij(len_gamma,:);
        iAtgxAgx_ij(len_gamma,:) = temp_row;
        temp_col = iAtgxAgx_ij(:,rowi);
        iAtgxAgx_ij(:,rowi) = iAtgxAgx_ij(:,len_gamma);
        iAtgxAgx_ij(:,len_gamma) = temp_col;
        
        AtgxAgx = AtgxAgx_ij(1:len_gamma-1,1:len_gamma-1);
        
        %iAtgxAgx = update_inverse(AtgxAgx_ij, iAtgxAgx_ij,2);
        n = size(AtgxAgx_ij,1);
        %delete columns
        Q11 = iAtgxAgx_ij(1:n-1,1:n-1);
        Q12 = iAtgxAgx_ij(1:n-1,n);
        Q21 = iAtgxAgx_ij(n,1:n-1);
        Q22 = iAtgxAgx_ij(n,n);
        Q12Q21_Q22 = Q12*(Q21/Q22);
        iAtgxAgx = Q11 - Q12Q21_Q22;
        
        xk_1(out_x(1)) = 0;
    else
        % If an element is added to gamma_x
        gamma_xk = [gamma_x; i_delta];
        new_x = i_delta;

        AtgxAnx = A(:,gamma_x)'*A(:,new_x);
        AtgxAgx_mod = [AtgxAgx AtgxAnx; AtgxAnx' A(:,new_x)'*A(:,i_delta)];
        
        AtgxAgx = AtgxAgx_mod;
        
        %iAtgxAgx = update_inverse(AtgxAgx, iAtgxAgx,1);
        n = size(AtgxAgx,1);
        % add columns
        iA11 = iAtgxAgx;
        iA11A12 = iA11*AtgxAgx(1:n-1,n);
        A21iA11 = AtgxAgx(n,1:n-1)*iA11;
        S = AtgxAgx(n,n)-AtgxAgx(n,1:n-1)*iA11A12;
        Q11_right = iA11A12*(A21iA11/S);
        %     Q11 = iA11+ Q11_right;
        %     Q12 = -iA11A12/S;
        %     Q21 = -A21iA11/S;
        %     Q22 = 1/S;
        iAtgxAgx = zeros(n);
        %iAtB = [Q11 Q12; Q21 Q22];
        iAtgxAgx(1:n-1,1:n-1) = iA11+ Q11_right;
        iAtgxAgx(1:n-1,n) = -iA11A12/S;
        iAtgxAgx(n,1:n-1) =  -A21iA11/S;
        iAtgxAgx(n,n) = 1/S;
        
        xk_1(i_delta) = 0;
    end

    z_xk = zeros(N,1);
    z_xk(gamma_xk) = -sign(Primal_constrk(gamma_xk));
    Primal_constrk(gamma_x) = sign(Primal_constrk(gamma_x))*epsilon;
end
total_iter = iter;

x_out = xk_1;

timeSteps = timeSteps(1:total_iter) ;
errorSteps = errorSteps(1:total_iter) ;
epsSteps = epsSteps(1:total_iter);

% Debiasing
%x_out = zeros(N,1);
%x_out(gamma_x) = A(:,gamma_x)\b;

% update_primal.m
%
% This function computes the minimum step size in the primal update direction and
% finds change in the primal or dual support with that step.
%
% Inputs:
% gamma_x - current support of x
% gamma_lambda - current support of lambda
% z_x - sign sequence of x
% z_lambda - sign sequence of lambda
% del_x_vec - primal update direction
% pk_temp
% dk
% epsilon - current value of epsilon
% out_lambda - element removed from support of lambda in previous step (if any)
%
% Outputs:
% i_delta - index corresponding to newly active primal constraint (new_lambda)
% out_x - element in x shrunk to zero
% delta - primal step size
%
% Written by: Salman Asif, Georgia Tech
% Email: sasif@ece.gatech.edu

function [i_delta, delta, out_x] = update_primal(out_x)

global N gamma_x z_x  xk_temp del_x_vec pk_temp dk epsilon isNonnegative

gamma_lc = setdiff(1:N, union(gamma_x, out_x));

if isNonnegative
    delta1 = inf;
else
    delta1_constr = (epsilon-pk_temp(gamma_lc))./(1+dk(gamma_lc));
    delta1_pos_ind = find(delta1_constr>0);
    delta1_pos = delta1_constr(delta1_pos_ind);
    [delta1 i_delta1] = min(delta1_pos);
    if isempty(delta1)
        delta1 = inf;
    end
end

delta2_constr = (epsilon+pk_temp(gamma_lc))./(1-dk(gamma_lc));
delta2_pos_ind = find(delta2_constr>0);
delta2_pos = delta2_constr(delta2_pos_ind);
[delta2 i_delta2] = min(delta2_pos);
if isempty(delta2)
    delta2 = inf;
end

if delta1>delta2
    delta = delta2;
    i_delta = gamma_lc(delta2_pos_ind(i_delta2));
else
    delta = delta1;
    i_delta = gamma_lc(delta1_pos_ind(i_delta1));
end

delta3_constr = (-xk_temp(gamma_x)./del_x_vec(gamma_x));
delta3_pos_index = find(delta3_constr>0);
[delta3 i_delta3] = min(delta3_constr(delta3_pos_index));
out_x_index = gamma_x(delta3_pos_index(i_delta3));

out_x = [];
if ~isempty(delta3) && (delta3 > 0) && (delta3 <= delta)
    delta = delta3;
    out_x = out_x_index;
end

%%% THESE ARE PROBABLY UNNECESSARY 
%%% NEED TO REMOVE THEM. 

% The following checks are just to deal with degenerate cases when more
% than one elements want to enter or leave the support at any step
% (e.g., Bernoulli matrix with small number of measurements)

% This one is ONLY for those indices which are zero. And we don't know where
% will its dx point in next steps, so after we calculate dx and its in opposite
% direction to z_x, we will have to remove that index from the support.
xk_1 = xk_temp+delta*del_x_vec;
xk_1(out_x) = 0;

wrong_sign = find(sign(xk_1(gamma_x)).*z_x(gamma_x)==-1);
if isNonnegative
    wrong_sign = union(wrong_sign, find(xk_1(gamma_x)<0));
end
if ~isempty(gamma_x(wrong_sign))
    delta = 0;
    % can also choose specific element which became non-zero first but all
    % that matters here is AtA(gx,gl) doesn't become singular.
    [val_wrong_x ind_wrong_x] =  sort(abs(del_x_vec(gamma_x(wrong_sign))),'descend');
    out_x = gamma_x(wrong_sign(ind_wrong_x));
end

% If more than one primal constraints became active in previous iteration i.e.,
% more than one elements wanted to enter the support and we added only one.
% So here we need to check if those remaining elements are still active.
i_delta_temp = gamma_lc(abs(pk_temp(gamma_lc)+delta*dk(gamma_lc))-(epsilon-delta) >= 10*eps);
if ~isempty(i_delta_temp)

    i_delta_more = i_delta_temp;
    if (length(i_delta_more)>=1) && (~any((i_delta_temp==i_delta)))
        % ideal way would be to check that incoming element doesn't make AtA
        % singular!
        [v_temp i_temp] = max(-pk_temp(i_delta_more)./dk(i_delta_more));
        i_delta = i_delta_more(i_temp);
        delta = 0;
        out_x = [];
    end
end