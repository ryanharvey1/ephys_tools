function [av, sd, sk, kt, cv] = waveStats(tt_tsd, norm)
%
% [av, sd, sk, kt, cv] = waveStats(tt_tsd, norm)
%
% calcualate basic waveform stats for each channel:
% 
% INPUT: 
%    tt_tsd  ... a TT file tsd with data array: Data(tt_tsd) = nSpikes x nCh x nSamp 
%    norm    ... 1 = normalize waveforms to |w|=1;  0 = don't normalize.     
%
% OUTPUT: 
%   av{nCh}(nSamp)  .... cell array of nSamp x 1 average waveforms array
%   sd{nCh}(nSamp)  .... cell array of nSamp x 1 standard deviations array 
%   sk{nCh}(nSamp)  .... cell array of nSamp x 1 skewness array 
%   kt{nCh}(nSamp)  .... cell array of nSamp x 1 kurtosis array 
%   cv{nCh}(nSamp x nSamp) 
%                   .... cell array of nSamp x nSamp waveform covariance matrices

% CCC: PL

ttd = Data(tt_tsd);
[nSpikes, nCh, nSamp] = size(ttd);

if norm
   % normalize waveforms to unit L2 norm (so that only their SHAPE or
   % relative angles but not their length (energy) matters)
   l2norms = sqrt(sum(ttd.^2,3));
   ttd = ttd./l2norms(:,:,ones(1,nSamp));
end%if

tmpsk = zeros(nSamp,1);
tmpkt = zeros(nSamp,1);

for iC = 1:nCh; 
   av{iC}=mean(squeeze(ttd(:,iC,:)))';
   cv{iC}=cov(squeeze(ttd(:,iC,:))); 
   sd{iC}=sqrt(diag(cv{iC}));
   for iS = 1:nSamp
      tmpsk(iS,1)=skewness(squeeze(ttd(:,iC,iS)));
      tmpkt(iS,1)=kurtosis(squeeze(ttd(:,iC,iS)));
   end%for
   sk{iC}=tmpsk;
   kt{iC}=tmpkt;
end%for