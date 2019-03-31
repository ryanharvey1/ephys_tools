function weights = ComputeXCPrediction(q,Q,ep,varargin)

% weights = ComputePeerPrediction(S,Q,ep)
% S: ts of the cell
% Q: binned spike train matrix (tsd object) of all other cells

q   = full(Data(Restrict(q,ep)));
Q   = Restrict(Q,ep);
rg  = Range(Q);
dt  = median(diff(rg));

pop = full(Data(Q)); %In Hz
if ~isempty(varargin)
    weights0 = varargin{1};
else
    weights0 = randn(size(pop,2),1);
end

%q_g     = gpuArray(q);
%pop_g   = gpuArray(pop);

%L = @(x)-SpkTrainLogLikelihood(q_g,dt*modifiedExp(pop_g*x)) + 0.25*x'*x;
L = @(x)-SpkTrainLogLikelihood(q,dt*modifiedGaussConv(pop*x));
options = optimoptions('fminunc','display','off');

problem.objective   = L;
problem.options     = options;
problem.solver      = 'fminunc';
problem.x0          = weights0;
problem.ObjectiveLimit = 1e-10;

weights = fminunc(problem);