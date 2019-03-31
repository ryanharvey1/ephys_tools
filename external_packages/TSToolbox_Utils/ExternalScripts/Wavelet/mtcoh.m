% MTCOH         multi-channel coherency (not time-resolved)
% 
% [ c, freq ] = mtcoh( varargin )
%
% based on mtcsd.m (see help there)
% returns coherency (complex) 
% when called without output arguments, plots spectra, coherence, and phase
%   of the last pair
%
% coherence = abs( c );
% phase = angle( c );
%
% calls: mtparam, mtcsd

% 04-nov-11 ES

function [ c, freq ] = mtcoh( varargin )

% Multitaper coherence density
%
% This basically does the same thing as mtcsd, but scales
% the cross-spectra to become coherences

%[x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers] = mtparam(varargin);

[y freq] = mtcsd(x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers);

freq = freq( fidx );
nCh1 = size(y,2);
nCh2 = size(y,3);

c = zeros(size(y));
for Ch1 = 1:nCh1
	for Ch2 = 1:nCh2
        c(:,Ch1, Ch2) = y(:,Ch1, Ch2) ./ sqrt( y (:, Ch1, Ch1 ) .* y(:, Ch2, Ch2 ) );
	end
end

% plot
if nargout < 1
    Ch1 = 1; 
    Ch2 = 2;
    Sxx = y( :, Ch1, Ch1 );
    Syy = y( :, Ch2, Ch2 );
    C = abs( c( :, Ch1, Ch2 ) );
    p = angle( c( :, Ch1, Ch2 ) );
    
    figure
    subplot(3,1,1);
    set(gca,'box','on');
    lx = line( freq, Sxx );
    ly = line( freq, Syy );
    set(lx,'color','k');
    set(ly,'color','r');
    %set(gca,'xscale','log','yscale','log');
    %set(gca,'xscale','linear','yscale','log');
    xlabel('frequency (Hz)');
    ylabel('power spectra (1/s)');
    legend('x','y');
    xlim( freq( [ 1 end ] ) );
    %    ylim([5,1e3]);
    
    subplot(3,1,2);
    set(gca,'box','on');
    l=line(freq,C);
    set(l,'color','k');
    %set(gca,'xscale','log','yscale','log');
    %set(gca,'xscale','linear','yscale','log');
    xlabel('frequency (Hz)');
    ylabel('coherence');
    xlim( freq( [ 1 end ] ) );
    
    subplot(3,1,3);
    set(gca,'box','on');
    l=line(freq,p);
    set(l,'color','k');
    %set(gca,'xscale','log','yscale','linear');
    xlabel('frequency (Hz)');
    ylabel('phase');
    xlim(freq( [ 1 end ] ) );
    ylim([-pi pi]);
end

return
