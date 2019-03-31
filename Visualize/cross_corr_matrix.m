function cross_corr_matrix(spikearray,spikesID)
% Plots CCG for all cell pairs
%
%   Input:
%           spikearray: cell array of spike times in seconds
%           spikesID: structure containing spikesID.TetrodeNum &
%                       spikesID.CellNum. Both the same length as your 
%                       spike array
%
%       
% Ryan Harvey 2018
fig = gcf;
rc=length(spikearray);
cor_store=[];
index = reshape(1:rc^2, rc, rc).';

for i=1:rc
    for ii=i:rc
        max_lag = 0.5;
        % t_bin=0.005;
        t_bin=0.01;
        % Acor - taken from intrinsic frequency 2
        if t_bin / mod(max_lag, t_bin) ~= 2 % set lags so it is 'even' (odd number of coefficients and zero centered')
            max_lag = t_bin*floor(max_lag/t_bin)+.5*t_bin;
        end
        [cor, lag] = CrossCorr(spikearray{i},spikearray{ii},'lag',[-max_lag max_lag],'binsize',t_bin,'norm','prob');
        
        % save cor
        cor_store=[cor_store,cor];
        
%         fig;
%         subplot(rc,rc,index(ii,i))
%         plot(lag,cor,'k')
%         box off
%         
%         grid off
%         title([erase(spikesID.TetrodeNum{i},'.mat'),'c',num2str(spikesID.CellNum(i)),' vs. ',...
%             erase(spikesID.TetrodeNum{ii},'.mat'),'c',num2str(spikesID.CellNum(ii))])
%         pause(.0001)
    end
end

cor_num=1;
for i=1:rc
    for ii=i:rc
        fig;
        subplot(rc,rc,index(ii,i))
        plot(lag,cor_store(:,cor_num),'k')
        cor_num=cor_num+1;
        box off
        
        grid off
        title([erase(spikesID.TetrodeNum{i},'.mat'),'c',num2str(spikesID.CellNum(i)),' vs. ',...
            erase(spikesID.TetrodeNum{ii},'.mat'),'c',num2str(spikesID.CellNum(ii))])
        pause(.0001)
    end
end