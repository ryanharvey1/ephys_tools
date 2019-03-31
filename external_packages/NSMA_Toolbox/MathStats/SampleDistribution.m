function S = SampleDistribution(bins, PDF, n)

% SampleDistribution  Samples a given distribution so one can do statistics on it
%
% S = SampleDistribution(bins, PDF, n)
% S = SampleDistribution(bins, PDF)    (n = 1000)
% S = SampleDistribution(PDF)          (bins = 1:length(PDF), n = 1000)
%
% INPUTS:
%       bins = vector representing binstarts for the PDF x-dimension.
%       PDF = vector representing counts for the PDF.
%       n = number of points to sample, default is 1000.
% OUTPUTS:
%       S = n points such that P(S in PDFx) = PDFy
%
% ADR 1998, version L4.0, last modified '98 by ADR

% status PROMOTED


switch (nargin)
case 1
   n = 1000;
   PDF = bins;
   bins = 1:length(PDF);
case 2
   n = 1000;
case 3
otherwise
   error('Incorrect number of arguments.');
end

   
sdf = sum(PDF);cdf = cumsum(PDF);

S = rand(n,1) * sdf;
for iT = 1:length(S)
   S(iT) = bins(binsearch(cdf, S(iT)));
end
