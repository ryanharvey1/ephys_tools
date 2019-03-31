% ThetaPower_testing
% THETA POWER INCREASE OVER LAPS?
% data=load('HPCatn02_S20180728180258')
% signal=mean(data.lfp.signal);
% for i=1:8
%     [MeanThetaFreq(i),MeanOverallPow(i)]=meanfreq(data.lfp.theta(i,:),1000);
% end
% [~,I]=max(MeanOverallPow);
I=3;
signal=data.lfp.signal(I,:);
theta=data.lfp.theta(I,:);  


for i=1:length(data.linear_track.right{1, 1}.laps)
    startr(i)=data.linear_track.right{1, 1}.laps{i}(1,1);
    endr(i)=data.linear_track.right{1, 1}.laps{i}(end,1);
end

for i=1:length(data.linear_track.left{1, 1}.laps)
    startl(i)=data.linear_track.left{1, 1}.laps{i}(1,1);
    endl(i)=data.linear_track.left{1, 1}.laps{i}(end,1);
end

for i=1:length(startr)
    lap_sig_r{i}=signal(data.lfp.ts>startr(i) & data.lfp.ts<endr(i));
end

for i=1:length(startl)
    lap_sig_l{i}=signal(data.lfp.ts>startl(i) & data.lfp.ts<endl(i));
end

[p,n]=numSubplots(length(lap_sig_r));
figure;
place=1;
for i=1:length(lap_sig_r)
    subplot(p(1),p(2),place)
    pspectrum(lap_sig_r{i},1000,'spectrogram','FrequencyLimits',[0, 60],'TimeResolution', 1);
    colormap jet
    place=place+1;
end

[p,n]=numSubplots(length(lap_sig_l));
figure;
place=1;
for i=1:length(lap_sig_l)
    subplot(p(1),p(2),place)
    pspectrum(lap_sig_l{i},1000,'spectrogram','FrequencyLimits',[0, 60],'TimeResolution', 1);
    colormap jet
    place=place+1;
end


%% mean theta per lap
for i=1:length(startr)
    [~,lap_theta_r(i)]=meanfreq(theta(data.lfp.ts>startr(i) & data.lfp.ts<endr(i)),1000);
end

for i=1:length(startl)
    [~,lap_theta_l(i)]=meanfreq(theta(data.lfp.ts>startl(i) & data.lfp.ts<endl(i)),1000);
end
figure;
subplot(1,2,2)
plot(1:length(lap_theta_r),lap_theta_r,'.k')
xlabel('trial')
ylabel('theta power')
lsline
[R,P] = corrcoef(1:length(lap_theta_r)',lap_theta_r) 


subplot(1,2,1)
plot(1:length(lap_theta_l),lap_theta_l,'.k')
xlabel('trial')
ylabel('theta power')
lsline
[R,P] = corrcoef(1:length(lap_theta_l)',lap_theta_l) 



% for i=1:length(lap_sig_r)
%     downsampledsigr(i,:)=interp1(1:length(lap_sig_r{i}),lap_sig_r{i},linspace(1,length(lap_sig_r{i}),1000));
% end
% mean(downsampledsigr);
% 
% 
% for i=1:length(lap_sig_l)
%     downsampledsigl(i,:)=interp1(1:length(lap_sig_l{i}),lap_sig_l{i},linspace(1,length(lap_sig_l{i}),1000));
% end
% mean(downsampledsigl);

% pspectrum(mean(data.lfp.signal),1000,'spectrogram','FrequencyLimits',[0, 60],'TimeResolution', 1);