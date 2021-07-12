function [ Z, E ] = low_rank( X, lambda, maxIter, A )
%LOW_RANK Summary of this function goes here
%   Detailed explanation goes here

d=size(X,1);
n=size(X,2);

Z=zeros(n,n);
J=Z;

E=zeros(d,n);
Y1=zeros(d,n);
Y2=zeros(n,n);

mu=1e-6;
max_mu=10^10;
rho=1.1;
epsilon=1e-8;

xtx=X'*X;
inv_x=inv(xtx+eye(n));

MAX_ITER=maxIter;
iter=0;

while true
    if iter>MAX_ITER
        fprintf(1,'max iter num reached!\n');
        h=figure('Visible', 'off');
        imagesc(Z);
        colormap(gray);
        axis equal;
        saveas(h,'tmp.png');
        break;
    end
    % 1. update J
    Y=Z+Y2/mu;
    tau=1/mu;
    J=singular_value_shrinkage(Y,tau);
    % 2. update Z
    Z=inv_x*(xtx-X'*E+J+(X'*Y1-Y2)/mu);
    Z=max(Z,0); % 非负矩阵
    Z=bsxfun(@times, Z, A); % pair-wise multiple, A is the spatial constraints matrix
    % 3. update E
    tlambda=lambda/mu;
    Q=X-X*Z+Y1/mu;
    E=l21(Q,tlambda);
    % 4. update Y1, Y2
    xz=X-X*Z;
    zj=Z-J;
    leq1=xz-E;
    leq2=zj;
    Y1=Y1+mu*(leq1);
    Y2=Y2+mu*(leq2);
    % 5. update mu
    mu=min(rho*mu,max_mu);
    % 6. check the convergence
    if max(max(abs(leq1)))<epsilon && max(max(abs(leq2)))<epsilon
        fprintf(1,'iter %d, convergenced\n', iter);
        h=figure('Visible', 'off');
        imagesc(Z);
        colormap(gray);
        axis equal;
        saveas(h,'tmp.png');
        % disp('press any key to continue');
        % pause;
        break;
    end
    if mod(iter, 50)==0
        % h=figure('Visible', 'off');
        % imagesc(Z);
        % colormap(gray);
        % axis equal;
        % saveas(h,'tmp.png');
        disp(['iter ' num2str(iter) ' press any key to continue']);
        % pause;
    end
    iter=iter+1;
end

end

