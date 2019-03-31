function [xs,Ws] = brownDecoder(Y,params,W)
%R is the response matrix (Ntimesteps x ncells)
%params are the parameters of the cells (ncells x 5, corresponding to
%x, y, sigmax, sigmay, and offset)
 
%W is the covariance matrix from p(x_k|x_k-1) = N(x_k-1,W^-1)
xhat = [0;0];
What = eye(2)*5;
 
xs = zeros(size(Y,1),2);
Ws = zeros(size(Y,1),2,2);
 
mfopts.Method = 'newton';
mfopts.Display = 'off';
mfopts.TolX = 1e-4;
mfopts.TolFun = 1e-4;
 
for ii = 1:size(Y,1)
%Computing p(x(t_k)|spike in (0,t_k-1)) is easy, it's simply
%N(xhat, What + W)
What = What + W;
 
%Next, figure out p(x_k|y_1...y_k) \approx N(xhat, What)
%through optimization
y = Y(ii,:)';
Wihat = inv(What);
xhat = minFunc(@(x) computeStatePosterior(x,y,params,xhat,Wihat),xhat,mfopts);
[~,~,Wihat] = computeStatePosterior(xhat,y,params,xhat,Wihat);
What = safeInverse(Wihat);
xs(ii,:) = xhat';
Ws(ii,:,:) = What;
 
if mod(ii,100) == 0
fprintf('Iteration %d\n',ii);
end
end
end
 
function [Hi] = safeInverse(H)
[~,p] = chol(H);
if p > 0
H = H + eye(size(H,1)) * max(0,10 - min(real(eig(H))));
end
Hi = inv(H);
end
 
function [E,g,H] = computeStatePosterior(x,y,params,xhat,Wihat)
%Predicted rate for each neuron given Gaussian RF models
loglambdas = -1/2./params(:,3).^2.*(x(1)-params(:,1)).^2 + ...
-1/2./params(:,4).^2.*(x(2)-params(:,2)).^2 + params(:,5);
lambdas = exp(loglambdas);
 
%Compute negative log-posterior
%First part is due to likelihood of data given position, second part is
%the prior prob. of positions
E = sum(-y.*loglambdas + lambdas) + 1/2*(x-xhat)'*Wihat*(x-xhat);
 
if nargout > 1
%Compute gradient of error (first part: likelihood)
dloglambdas = [-1./params(:,3).^2.*(x(1)-params(:,1)),-1./params(:,4).^2.*(x(2)-params(:,2))];
g1 = zeros(2,1);
g1(1) = sum(-y.*dloglambdas(:,1) + lambdas.*dloglambdas(:,1));
g1(2) = sum(-y.*dloglambdas(:,2) + lambdas.*dloglambdas(:,2));
%Second part: prior
g2 = Wihat*(x-xhat);
 
g = g1+g2;
end
 
if nargout > 2
 
%Compute Hessian of error, first part
H1(1,1) = sum(-y.*(-1./params(:,3).^2) + lambdas.*dloglambdas(:,1).^2 + lambdas.*(-1./params(:,3).^2));
H1(1,2) = sum(                          lambdas.*dloglambdas(:,1).*dloglambdas(:,2));
H1(2,1) = H1(1,2);
H1(2,2) = sum(-y.*(-1./params(:,4).^2) + lambdas.*dloglambdas(:,2).^2 + lambdas.*(-1./params(:,4).^2));
 
H2 = Wihat;
H = H1 + H2;
 
end
end