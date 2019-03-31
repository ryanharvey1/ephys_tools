function [R, decordist, oppcorr]=PopVecCorr_bidir(Q,plotflag,bins,titlestr)
% create correlation matrix from Q matrix (of population vectors)
% [R, decordist, oppcorr]=PopVecCorr_bidir(Q,plotflag,bins,titlestr)
% inputs: Q: the Q matrix, must have cells in rows, and time (or position)
%            bins in columns
%         bins: vector of the value of each time bin. for ploting of
%               correlation matrix
%         titlestr: description of Qmtx (title for plot)
%         plotflag: plot the correlation matrix?
% outputs: R: the population vector correlation matrix
%          decordist: average correlation along diagonal: the profile for
%                     bins to decorrelate
%          oppcorr: correlation of running in opposite direction
% ZN 02/2010
% last edited: 04/11

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
    title(titlestr)
    axis square
end

% find decorrelation distance
decordist=zeros(nBins, 1);
oppcorr=zeros((nBins*2)-1,1);
for a=0:nBins-1
    decordist(a+1)=nanmean(diag(R, a));
    oppcorr(a+nBins)=nanmean(diag(fliplr(R),a));
    oppcorr(nBins-a)=nanmean(diag(fliplr(R),-1*a));
end

if plotflag
    binsize=bins(2)-bins(1);
    figure; plot(-1*binsize*(nBins-1):binsize:binsize*(nBins-1), [flipud(decordist(2:end));decordist], 'b')
    hold on; plot(-1*binsize*(nBins-1):binsize:binsize*(nBins-1), oppcorr, 'r')
    xlabel('Distance')
    ylabel('Correlation')
    legend('normal', 'reverse direction')
    title(['Decorrelation distance: ', titlestr])
    axis tight
end


