function [P,sumP] = directIntersect(X,Y)
%UNTITLED2 Summary of this function goes here
%   depends on InterX
%path
L1=[X'; Y'];

%path same size as observed from start to end
xDirect = linspace(X(1,1),X(end,1),size(X,1)); % create X 
yDirect = linspace(Y(1,1),Y(end,1),size(X,1)); % create Y
L2=[xDirect; yDirect]; %concatenate for Interx

P = InterX(L1(:,2:(end-1)),L2(:,2:(end-1))); %index of intersections
sumP=size(P,2);

end
