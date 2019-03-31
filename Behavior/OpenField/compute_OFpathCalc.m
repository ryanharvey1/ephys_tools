function [ PathLength, pathIV] = compute_OFpathCalc(x_smooth,y_smooth,fr)
%compute_pathCalc calulates the length of a path given 
%   Detailed explanation goes here
x_smooth(isnan(x_smooth))=0;
y_smooth(isnan(y_smooth))=0;
sqrXDiff=(diff(x_smooth)).^2;
sqrYDiff=(diff(y_smooth)).^2;
pathDist=sqrt(sqrXDiff + sqrYDiff);                                         %distance formula
pathIV=(pathDist)*fr; % .10 = a frame rate of 10 hz                         %instanteous velocity
PathLength=sum(pathDist);                                                   %Total Path length

end

