% % phaseprecessiontesting
% load('/Users/ryanharvey/Downloads/phaseprec_workspace.mat')
% %% example phase 
% % https://www.gaussianwaves.com/2017/04/extracting-instantaneous-amplitude-phase-frequency-hilbert-transform/
% %
% fs = 600; %sampling frequency in Hz
% t = 0:1/fs:1-1/fs; %time base
% a_t = 1.0 + 0.7 * sin(2.0*pi*3.0*t) ; %information signal
% c_t = chirp(t,20,t(end),80); %chirp carrier
% x = a_t .* c_t; %modulated signal
%  
% subplot(2,1,1); plot(x);hold on; %plot the modulated signal
%  
% z = hilbert(x); %form the analytical signal
% inst_amplitude = abs(z); %envelope extraction
% inst_phase = unwrap(angle(z));%inst phase
% inst_freq = diff(inst_phase)/(2*pi)*fs;%inst frequency
%  
% %Regenerate the carrier from the instantaneous phase
% regenerated_carrier = cos(inst_phase);
%  
% plot(inst_amplitude,'r'); %overlay the extracted envelope
% title('Modulated signal and extracted envelope'); xlabel('n'); ylabel('x(t) and |z(t)|');
% subplot(2,1,2); plot(cos(inst_phase));
% title('Extracted carrier or TFS'); xlabel('n'); ylabel('cos[\omega(t)]');
% %%
% fs = 1000; %sampling frequency in Hz
% % t = 0:1/fs:1-1/fs; %time base
% % a_t = 1.0 + 0.7 * sin(2.0*pi*3.0*t) ; %information signal
% % c_t = chirp(t,20,t(end),80); %chirp carrier
% % x = a_t .* c_t; %modulated signal
%  
% subplot(2,1,1); plot(EEGthetaData);hold on; %plot the modulated signal
%  
% z = hilbert(EEGthetaData); %form the analytical signal
% inst_amplitude = abs(z); %envelope extraction
% inst_phase = unwrap(angle(z));%inst phase
% inst_freq = diff(inst_phase)/(2*pi)*fs;%inst frequency
%  
% %Regenerate the carrier from the instantaneous phase
% regenerated_carrier = cos(inst_phase);
%  
% spikephase=interp1(EEG_DownSampledTimestamps,regenerated_carrier,data_video_spk(data_video_spk(:,6)==1,1));
% 
% plot(inst_amplitude,'r'); %overlay the extracted envelope
% title('Modulated signal and extracted envelope'); xlabel('n'); ylabel('x(t) and |z(t)|');
% subplot(2,1,2); plot(cos(inst_phase));
% title('Extracted carrier or TFS'); xlabel('n'); ylabel('cos[\omega(t)]');
% 
% %%
% figure;
% subplot(4,1,1)
% plot(EEG_DownSampledData(1:1000),'k')
% subplot(4,1,2)
% plot(EEGthetaData(1:1000),'k')
% subplot(4,1,3)
% plot(inst_phase(1:1000),'k')
% subplot(4,1,4)
% plot(regenerated_carrier(1:1000),'k')
% %%
% close all
% figure;
% subplot(2,1,1)
% plot(data_video_spk(:,2),data_video_spk(:,1),'k')
% hold on;
% scatter(data_video_spk(data_video_spk(:,6)==1,2),data_video_spk(data_video_spk(:,6)==1,1), 35, 'filled', 'r')
% subplot(2,1,2)
% h=scatter(data_video_spk(data_video_spk(:,6)==1,2),spikephase,35, 'filled', 'r')
% h.CData=abs(spikeangle); colormap jet
% xlim([0 500])
% %%
% figure;
% subplot(2,1,1)
% plot(data_video_spk(:,2),data_video_spk(:,1),'k')
% hold on;
% h=scatter(data_video_spk(data_video_spk(:,6)==1,2),data_video_spk(data_video_spk(:,6)==1,1), 35, 'filled', 'r');
% h.CData=spikeangle; colormap jet
% 
% subplot(2,1,2)
% 
% spikeangle=(rescale(spikephase,-180,180));
% h=scatter(data_video_spk(data_video_spk(:,6)==1,2),spikeangle, 'filled', 'r');
% xlim([0 500])
% ylim([-180 180])
% 
% RSquared=[];
% for iPH=0:359
%     mdl=fitlm(data_video_spk(data_video_spk(:,6)==1,2),circshift(spikeangle,iPH));
%     RSquared=[RSquared;mdl.Rsquared.Ordinary];
%     h.YData=circshift(spikeangle,iPH);
%     refreshdata
%     lsline
%     pause(.000001)
% end
% %%
% [circ_lin_corr pval slope_deg phi0_deg RR] = cl_corr(data_video_spk(data_video_spk(:,6)==1,2), rescale(spikephase,-180,180), 0, 1)
% 
% 
% 
% %%  spikes on phase moving plot
% % figure;
% % scatter(EEG_DownSampledTimestamps,ones(length(EEG_DownSampledTimestamps),1)+1)
% % hold on
% % scatter(data_video_spk(data_video_spk(:,6)==1,1),ones(length(data_video_spk(data_video_spk(:,6)==1,1)),1))
% 
% figure;
% x1=EEG_DownSampledTimestamps(1);
% x2=EEG_DownSampledTimestamps(1000);
% 
% 
% plot(EEG_DownSampledTimestamps,rescale(regenerated_carrier,-180,180),'k');hold on
% h=scatter(data_video_spk(data_video_spk(:,6)==1,1),spikeangle, 35, 'filled', 'r')
% h.CData=spikeangle;
% xlim([x1 x2])
% 
% for i=1:10:length(EEG_DownSampledTimestamps)
%     xlim([EEG_DownSampledTimestamps(1+i) EEG_DownSampledTimestamps(1000+i)])
%     refreshdata
%     pause(.0000001)
% end
% 
% 
% %% try velocity filtering more? 
% EEGthetaData = eegfilt(EEG_DownSampledData,NewSFreq,6,10);  % THETA
% 
% %%
% 	phase= mod(angle(hilbert(EEGthetaData)),2*pi);
%     spikephase=interp1(EEG_DownSampledTimestamps,phase,data_video_spk(data_video_spk(:,6)==1,1));
%     spikephase=rad2deg(spikephase);
% figure;
% subplot(2,1,1)
% plot(data_video_spk(:,2),data_video_spk(:,1),'k')
% hold on;
% scatter(data_video_spk(data_video_spk(:,6)==1,2),data_video_spk(data_video_spk(:,6)==1,1), 35, 'filled', 'r')
% subplot(2,1,2)
% h=scatter([data_video_spk(data_video_spk(:,6)==1,2)],spikephase)
% % h.CData=abs(spikephase); colormap jet
% xlim([0 500])
% %%
% 
%     
% % spkval=data_video_spk(data_video_spk(:,6)==1,:);
% % newdata_video_spk=data_video_spk(data_video_spk(:,5)>5,:);
% % 
% % spkval=spkval(spkval(:,5)>5,:);
% % 
% % spikephase=interp1(EEG_DownSampledTimestamps,regenerated_carrier,spkval(:,1));
% % 
% % close all
% % figure;
% % subplot(2,1,1)
% % plot(newdata_video_spk(:,2),newdata_video_spk(:,1),'k')
% % hold on;
% % scatter(spkval(:,2),spkval(:,1), 35, 'filled', 'r')
% % subplot(2,1,2)
% % h=scatter(spkval(:,2),spikephase,35, 'filled', 'r')
% % h.CData=abs(spikephase); colormap jet
% % xlim([0 500])
% 
% 
% 
% 
% % z = hilbert(EEGthetaData); %form the analytical signal
% % inst_amplitude = abs(z); %envelope extraction
% % inst_phase = unwrap(angle(z));%inst phase
% % inst_freq = diff(inst_phase)/(2*pi)*fs;%inst frequency
% %  
% % %Regenerate the carrier from the instantaneous phase
% % regenerated_carrier = cos(inst_phase);
% % 
% % spikephase=interp1(EEG_DownSampledTimestamps,regenerated_carrier,spkval(:,1));
% % 
% % 
% % figure;
% % subplot(2,1,1)
% % plot(newdata_video_spk(:,2),newdata_video_spk(:,1),'k')
% % hold on;
% % scatter(spkval(:,2),spkval(:,1), 35, 'filled', 'r')
% % subplot(2,1,2)
% % h=scatter(spkval(:,2),regenerated_carrier,35, 'filled', 'r')
% % h.CData=abs(regenerated_carrier); colormap jet
% % xlim([0 500])
% 
% 
% % h=scatter(spkval(1:5,2),spikephase(1:5),35, 'filled', 'r')
% % figure;comet(spkval(1:5,2),spikephase(1:5))
% 
% % % make it move
% % figure;
% % 
% % subplot(2,1,1)
% % h1=plot(NaN,'k');
% % hold on
% % h2=scatter(NaN,NaN, 35, 'filled', 'r');
% % xlim([0 500])
% % ylim([min(newdata_video_spk(:,1)) max(newdata_video_spk(:,1))])
% % 
% % subplot(2,1,2)
% % h3=scatter(NaN,NaN, 35, 'filled', 'k');
% % xlim([0 500])
% % ylim([-1 1])
% % 
% % for i=1:length(newdata_video_spk(:,1))
% %     h1.XData=newdata_video_spk([1:i],2);
% %     h1.YData=newdata_video_spk([1:i],1);
% %     if newdata_video_spk(i,6)==1
% %         
% %         h2.XData=spkval([1:i],2);
% %         h2.YData=spkval([1:i],1);
% %         
% %         h3.XData=spkval([1:i],2);
% %         h3.YData=spikephase(1:i);
% %         refreshdata
% %         pause(1)
% %     end
% %     
% % end
% 
% %% new test 1/30/18
% 
% fs = 1000; %sampling frequency in Hz
% % t = 0:1/fs:1-1/fs; %time base
% % a_t = 1.0 + 0.7 * sin(2.0*pi*3.0*t) ; %information signal
% % c_t = chirp(t,20,t(end),80); %chirp carrier
% % x = a_t .* c_t; %modulated signal
%  
% subplot(2,1,1); plot(EEGthetaData);hold on; %plot the modulated signal
%  
% z = hilbert(EEGthetaData); %form the analytical signal
% inst_amplitude = abs(z); %envelope extraction
% inst_phase = unwrap(angle(z));%inst phase
% inst_freq = diff(inst_phase)/(2*pi)*fs;%inst frequency
%  
% %Regenerate the carrier from the instantaneous phase
% regenerated_carrier = cos(inst_phase);
%  
% spikephase=interp1(EEG_DownSampledTimestamps,regenerated_carrier,data_video_spk(data_video_spk(:,6)==1,1));
% spikephase2=spline(EEG_DownSampledTimestamps,regenerated_carrier,data_video_spk(data_video_spk(:,6)==1,1));
% 
% % max
% % min
% rescale(regenerated_carrier,0,2*pi);
% 
% % phiOut = wrapToPi(cumsum(abs(diff(regenerated_carrier))));
% 
% 
% %  exp(i*regenerated_carrier);
% spikephase=angle(interp1(EEG_DownSampledTimestamps,exp(1i*rescale(regenerated_carrier,0,2*pi)),data_video_spk(data_video_spk(:,6)==1,1)));
% angle_rad = spikephase - 2*pi*floor(spikephase/(2*pi));
% %%
% ThetaPh = angle(hilbert(EEGthetaData));
% 
% spikephase=interp1(EEG_DownSampledTimestamps,ThetaPh,data_video_spk(data_video_spk(:,6)==1,1));
% 
% spkval=data_video_spk(data_video_spk(:,6)==1,:);
% 
% spikesums=cumsum(data_video_spk(:,6));
% %%
% for i=1:length(spkval)
% nspks=i;
% figure;
% subplot(2,1,1)
% plot(data_video_spk([1:find(spikesums==nspks)],2),data_video_spk([1:find(spikesums==nspks)],1))
% hold on
% scatter(spkval([1:nspks],2),spkval([1:nspks],1),'r')
% xlim([0 500])
% 
% subplot(2,1,2)
% h=scatter(spkval([1:nspks],2),spikephase(1:nspks),35, 'filled', 'r')
% xlim([0 500])
% pause(2)
% close all
% end
%% 1/31/18
% ThetaPh = angle(hilbert(EEGthetaData));
% spikephase=interp1(EEG_DownSampledTimestamps,ThetaPh,data_video_spk(data_video_spk(:,6)==1,1));
% figure;scatter(data_video_spk(data_video_spk(:,6)==1,2),spikephase)

% spikephase=interp1(EEG_DownSampledTimestamps,EEGthetaData,data_video_spk(data_video_spk(:,6)==1,1));
%%
% figure;scatter(data_video_spk(data_video_spk(:,6)==1,2),spikephase,'k')
% hold on
% scatter(data_video_spk(data_video_spk(:,6)==1,2),spikephase+max(spikephase),'k')

% ThetaPh = angle(hilbert(EEGthetaData));


% inst_phase = cos(unwrap(angle(hilbert(EEGthetaData))));

fs=1000;
z = hilbert(EEGthetaData); %form the analytical signal
inst_amplitude = abs(z); %envelope extraction
inst_phase = unwrap(angle(z));%inst phase
inst_freq = diff(inst_phase)/(2*pi)*fs;%inst frequency
%Regenerate the carrier from the instantaneous phase
regenerated_carrier = cos(inst_phase);

spikephase=interp1(EEG_DownSampledTimestamps,regenerated_carrier,data_video_spk(data_video_spk(:,6)==1,1));

figure;
h=scatter(data_video_spk(data_video_spk(:,6)==1,2),spikephase, 'filled', 'r');
xlim([0 500])

RSquared=[];
for iPH=0:359
    mdl=fitlm(data_video_spk(data_video_spk(:,6)==1,2),circshift(spikephase,iPH));
    RSquared=[RSquared;mdl.Rsquared.Ordinary];
    h.YData=circshift(spikephase,iPH);
    refreshdata
    lsline
    pause(.000001)
end

figure;plot(EEG_DownSampledTimestamps,regenerated_carrier)
hold on
h=scatter(data_video_spk(data_video_spk(:,6)==1,1),spikephase, 'filled', 'r')
%% 2/1/18
% 
% 
% function [ph, thpeaks] = ThetaPhase(S, CRtheta, ThStart, ThEnd)
% ph = ThetaPhase(S, CRtheta, ThStart, ThEnd)
%
% Computes the theta phase of each spike for a group of cells 
%
% INPUTS:
% S:              cell array of ts objects containing the spike trains
% CRtheta:        EEG tsd object filtered for theta
% ThStart, Thend: arrays containig the start and stop times of the valid
%                 theta epochs
%
% OUTPUTS: 
% ph:             cell array of tsd abjects containing the  timestamp of each
%                 spike in Range(ph,'ts') and the phase S (in range [0,1] and the theta cycle number
%                 in a 2 column array in Data(ph)
% thpeaks:        The timestamps of each theta peak (the peak sample point of CRtheta tsd)
%
% batta 2000 (modified by PL Aug. 2000)
% status: alpha
addpath(genpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/eeglab14_0_0b'))
addpath(genpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/NSMA_Toolbox'))

fs=1000;
EEGthetaData = eegfilt(EEG_DownSampledData,fs,4,12);  % THETA

% t = Range(CRtheta, 'ts');
t=EEG_DownSampledTimestamps';
% dd = Data(CRtheta);
dd=EEGthetaData';

dth = diff(dd);

dth1 = [0 dth'];
dth2 = [dth' 0];
clear dth;


peakindex = find(dth1 > 0 & dth2 < 0);
thpeaks = t(peakindex);
thpeaksAmplitudes = dd(peakindex);


clear t;
% ph = cell(length(S),1);
% for iC = 1:length(S)
%   SiC = Restrict(S{iC}, ThStart, ThEnd); 
%   s = Data(SiC);
%   s = Data(S{iC});
filtspks=spks_VEL(spks_VEL(:,5)>10,:);
  s=filtspks(:,1);
  
  
  ph = zeros(size(s));
  pks = zeros(size(s));
  for j = 1:length(s)
    pk = binsearch_floor(thpeaks, s(j));
    if pk ~= length(thpeaks)
    ph(j) = (s(j) - thpeaks(pk)) / (thpeaks(pk+1) - thpeaks(pk));
    pks(j) = pk;
    end
  end
%   ph = tsd( s, [ph pks]);
% end
% end
% [ph',pks'];
%
%  PLOT RESULTS FROM THE ABOVE METHOD
close all
figure;
subplot(2,1,1)
plot(data_video_spk(:,2),data_video_spk(:,1),'k')
hold on;
h=scatter(filtspks(:,2),filtspks(:,1), 35, 'filled', 'r');
% h.CData=abs(ph'); colormap jet

subplot(2,1,2)
plot(filtspks(:,2),ph','k.')
hold on
plot(filtspks(:,2),ph'+1,'b.')
xlim([0 500])
ylabel('Theta Phase')
xlabel('Position on Track')

% bin and plot ratemap
nBinsx = 50;
nBinsy = 50;
%     mat4occ=unique (occMatrix,'rows');
    MinY = min(ph');
    MaxY = max(ph'+1);
    MinX = 0;
    MaxX = 500;
    edges{1} = linspace(MinY, MaxY, nBinsy+1);
    edges{2} = linspace(MinX, MaxX, nBinsx+1);
    
    Omatrix = hist3([[ph;ph+1], [filtspks(:,2);filtspks(:,2)]],'Edges',edges);
    
    filtWidth = [5 5]; filtSigma = 1;
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    Omatrix = nanconv(Omatrix,imageFilter, 'nanout');
    
    figure;pcolor(Omatrix);shading interp;colormap jet
    axis xy
    
