function fixedBiasMatrix = BiasRemoveNaNs(biasMatrix)
%
% remove all rows in biasMatrix = [biasS1, biasB, biasS2] which have one or more NaNs in any of the S1, B or S2 
% epochs.
% returns a matrix fixedBiasMatrix with no NaNs.
%
% PL April 2003

ss = sum(isnan(biasMatrix),2);
ix_rows = find(ss);
fixedBiasMatrix = biasMatrix;
fixedBiasMatrix(ix_rows,:) = [];

