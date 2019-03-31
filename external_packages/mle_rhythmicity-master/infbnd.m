function in = infbnd(in)
%INFBND Handles bad parameter sets
%   Helper function that prevents infitite/negative likelihoods when 
%   parameters are bad.
    in(isinf(in)|isnan(in)|~isreal(in)) = -realmax;
end