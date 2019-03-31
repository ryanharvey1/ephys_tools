function [xs,Ws,particles] = pfDecoder(Y,params,W)
% https://xcorr.net/2012/02/07/using-a-particle-filter-to-decode-place-cells/

%Here I implement the method of Brockwell, Rojas and Kass (2004)
nparticles = 1000;
 
%Draw from the initial state distribution
xtildes = randn(nparticles,2);
 
xs = zeros(size(Y,1),2);
Ws = zeros(size(Y,1),2,2);
particles = zeros(size(Y,1),nparticles,2);
 
Wsqrt = W^(1/2);
 
%Main loop
for ii = 1:size(Y,1)
%Step 2: compute weights according to w = p(y|v_t=xtilde)
ws = computeLogLikelihoods(xtildes,Y(ii,:)',params);
%Normalize weights
ws = exp(ws-max(ws));
ws = ws/sum(ws);
 
%Step 3: importance sampling
idx = mnrnd(nparticles,ws);
S = bsxfun(@le,bsxfun(@plus,zeros(nparticles,1),1:max(idx)),idx');
[idxs,~] = find(S);
xtildes = xtildes(idxs,:);
 
particles(ii,:,:) = xtildes;
 
%Step 3.5
xs(ii,:) = mean(xtildes);
Ws(ii,:,:) = cov(xtildes);
 
%Step 4: Propagate each particle through the state equation
xtildes = xtildes + (Wsqrt*randn(2,nparticles))';
 
if mod(ii,100) == 0
fprintf('Iteration %d\n',ii);
end
 
end
end
 
function lls = computeLogLikelihoods(xtildes,y,params)
%Predicted rate for each neuron given Gaussian RF models
xdelta = bsxfun(@times,bsxfun(@minus,xtildes(:,1),params(:,1)'),1./(params(:,3)')).^2;
ydelta = bsxfun(@times,bsxfun(@minus,xtildes(:,2),params(:,2)'),1./(params(:,4)')).^2;
 
loglambdas = bsxfun(@plus,-.5*(xdelta+ydelta),params(:,5)');
lambdas = exp(loglambdas);
 
%Compute negative log-posterior
%First part is due to likelihood of data given position, second part is
%the prior prob. of positions
lls = -(-loglambdas*y + sum(lambdas,2));
end