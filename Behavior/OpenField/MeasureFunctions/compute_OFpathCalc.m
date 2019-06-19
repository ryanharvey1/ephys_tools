function [ PathLength, pathIV] = compute_OFpathCalc(x_smooth,y_smooth,fr)
%compute_pathCalc calulates the length of a path given 
%   Detailed explanation goes here
x_smooth(isnan(x_smooth))=0;
y_smooth(isnan(y_smooth))=0;
sqrXDiff=(diff(x_smooth)).^2;
sqrYDiff=(diff(y_smooth)).^2;
pathDist=sqrt(sqrXDiff + sqrYDiff);%distance formula
IV=(pathDist)*fr; % 10 = a frame rate of 10 hz      %instanteous velocity
pathIV=smoothdata(IV,'movmedian',fr*.8);

PathLength=sum(pathDist(pathIV>=3,1));%Total Path length for points greater than 3cm/s

end

