function [ccca, bincenters] = MakeCrossCorrCellArray (S, bnsz, nbins, sameTT)

%
% INPUTS:
%      S =      the "S-Matrix"; a cell array of ts-objects of spiketimes; one ts-object per cell.
%      bnsz =   the binsize for the cross corr histogram in msec.
%      nbins =  number of bins (timelags) in cross-correlograms;
%                   binsize * nbins+1 = total cross-correlograms time range in msec
%      sameTT = the same-NTrode-Matrix is a N*N matrix with ones for cellpairs from different probes and
%                   zeros for cellpairs from same probe; (e.g. the matrix returned by sameNtMatrix = getSameNTrodeMatrix(ses))
%
% OUTPUTS:
%      ccca = cross-correlogram cell array ccca{1..npairs} of  cross-correlograms (nbins+1 vectors);
%             only the lower triangle j>i of the symmetric matrix of cell-pairs (i,j) array is computed.
%             cell pairs with sameTT=1 are not included in the pairs-list;
%             The normalization of the the cross-corr histograms is different (multiplied by number of spikes in first cell)
%             from the CrossCorr.dll function!
%      bincenters = vector of bincenter coordinates of the cross-correlograms (which all have identical x-axes)
% 
% PL  April 18, 2003
%%%%%%%%%%%%%%%


NCells = length(S);
NPairs = 0;
ccca = {};
for ii = 1:NCells;
    for jj = ii+1:NCells;
        if isempty(Data(S{ii})) | isempty(Data(S{jj}))
            % don't count pairs where at least one cell is empty
            continue;
        end
        if sameTT(ii,jj) == 0
            % don't count pairs of cells on the same NTrode
            continue;
        end
        Na = length(Data(S{ii}));
        [tmp,bincenters] = CrossCorr(Data(S{ii}), Data(S{jj}),bnsz,nbins);
        tmp = tmp*Na;     % multiply by number of spikes of cell A to compensate for our
                          % normalization of CrossCorr, which implicitly divides by N1 to get
                          % the conditional firing rate of B;
                          % this way the CrossCorr between cell A and B is equal to the 
                          % flipped CrossCorr between cell B and A (PL).
         NPairs = NPairs+1;
         ccca{NPairs} = tmp;
    end %for jj
end %for ii
        
