function [R, decordist]=PopVecCorr(Q,plotflag,bins,titlestr)
% create correlation matrix from Q matrix
% [R, decordist]=PopVecCorr(Q,plotflag,bins,titlestr)
%
% inputs: Q: the Q matrix, must have cells in rows, and time (or position)
%            bins in columns
%         plotflag(optional): 1 if want plots, 0 if just want R and
%                             decordist (default is 1)
%         bins(optional): vector of the value of each time bin. for ploting of
%                         correlation matrix
%         titlestr(optional): description of Qmtx (title for plot)
% 
% outputs: R: the population vector correlation matrix
%          decordist: average correlation along diagonal: the profile for
%                     bins to decorrelate
% ZN 02/2010

[nCells,nBins]=size(Q);
if plotflag || nargin<2
    plotflag=1;
    if nargin<4
        titlestr='';
        if nargin<3
            bins=1:nBins;
        end
    end
    if length(bins)~=nBins
        disp('Warning: bins entered do not match bins in Q')
        bins=1:nBins;
    end
end

% autocorr the Q matrix
R=corrcoef_AB(Q, Q);

if plotflag
    figure; imagesc(bins, bins, R)
    title(['Population vector cross correlation: ', titlestr])
    axis square
end

% find decorrelation distance
decordist=zeros(nBins, 1);
for a=0:nBins-1
    decordist(a+1)=nanmean(diag(R, a));
end

if plotflag
    binsize=bins(2)-bins(1);
    figure; plot(0:binsize:binsize*(nBins-1), decordist, 'r')
    xlabel('Distance')
    ylabel('Correlation')
    title(['Decorrelation distance: ', titlestr])
    axis tight
end

