function [ a,b,c ] = min_encl_ellipsoid( X,Y )
%Finds the minimum enclosing ellipsoid of points along a rat's trajectory
%in the Morris Water Maze
%   X=all x-ccordinates of the trajectory
%   Y=all y-coordinates of the trajectory
%   Coordinates must be in columns, not rows
%   a=major radius of the ellipse
%   b=minor radius of the ellipse
%   c=center of the ellipse
% Finds the minimum volume enclsing ellipsoid (MVEE) of a set of data
% points stored in matrix P. The following optimization problem is solved: 
%
% minimize       log(det(A))
% subject to     (P_i - c)' * A * (P_i - c) <= 1
%                
% in variables A and c, where P_i is the i-th column of the matrix P. 
% The solver is based on Khachiyan Algorithm, and the final solution 
% is different from the optimal value by the pre-spesified amount of 'tolerance'.
%
%
% Nima Moshtagh (nima@seas.upenn.edu)
% University of Pennsylvania
%
% December 2005
% UPDATE: Jan 2009



%%%%%%%%%%%%%%%%%%%%% Solving the Dual problem%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% ---------------------------------
% data points 
% -----------------------------------
L=horzcat(X,Y);
P=L';
[d, N] = size(P);

Q = zeros(d+1,N);
Q(1:d,:) = P(1:d,1:N);
Q(d+1,:) = ones(1,N);


% initializations
% -----------------------------------
count = 1;
err = 1;
tolerance=.01;
u = (1/N) * ones(N,1);          % 1st iteration


% Khachiyan Algorithm
% -----------------------------------
while err > tolerance,
    X = Q * diag(u) * Q';       % X = \sum_i ( u_i * q_i * q_i')  is a (d+1)x(d+1) matrix
    M = diag(Q' * inv(X) * Q);  % M the diagonal vector of an NxN matrix
    [maximum, j] = max(M);
    step_size = (maximum - d -1)/((d+1)*(maximum-1));
    new_u = (1 - step_size)*u ;
    new_u(j) = new_u(j) + step_size;
    count = count + 1;
    err = norm(new_u - u);
    u = new_u;
end



%%%%%%%%%%%%%%%%%%% Computing the Ellipse parameters%%%%%%%%%%%%%%%%%%%%%%
% Finds the ellipse equation in the 'center form': 
% (x-c)' * A * (x-c) = 1
% It computes a dxd matrix 'A' and a d dimensional vector 'c' as the center
% of the ellipse. 

U = diag(u);

% the A matrix for the ellipse
% --------------------------------------------
A = (1/d) * inv(P * U * P' - (P * u)*(P*u)' );


% center of the ellipse 
% --------------------------------------------
c = P * u;
[~, S, ~,] = svd(A);
r1=1/sqrt(S(1,1));
r2=1/sqrt(S(2,2));
if r1>r2
    a=r1; b=r2;
else 
    a=r2; b=r1;
end
end

