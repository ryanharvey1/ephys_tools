%component_mu calculates the component means used in a von Mises Gaussian
%Mixture model analysis.
%
% created by Teagan Mullins Aug 2019
%
%Input:
%   numDist: number (double) indicating how many sets of means to create.

%Output:
%   circlepts: cell array (1 x numDist) containing nested (1 by numdist,
%   double) array of means in radians. Means are created to be equidistant given the
%    dist n (e.g. n=3 will have mu's of 0,120,240).
%

function circlepts=component_mu(numDist)


for i = 1:numDist
    spacing = 360/i;
    
    for j = 1:i
        n (j,1) = wrapToPi(deg2rad((j-1)* spacing));
    end
    circlepts{1,i} = n;
    
end

end