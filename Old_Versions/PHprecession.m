function [ ThPrecess ] = PHprecession(phase,spks_VEL,occ4Ph,fieldbound)
%PHprecession filters for theta and calculates phase precession
%
% Input:    EEG_DownSampledData:        raw downsampled lfp data
%           EEG_DownSampledTimestamps:  timestamps associated with downsampled lfp data
%           spks_VEL:                   timestamps associated with spike occurrences
%           NewSFreq:                   downsampled frequency
%           track_length:               length or diameter of apparatus
% ----------------------------------------------------------------------------------------
% Output:   ThPrecess
%               phaselock:
%                           Rlength:R length of spike phases
%                           Pval:   p value of R length
%               slope:              slope of regression line
%               RSquared:           R-Squared
%               scatteredPH:        phase by position matrix for scatter plot
%               lapSlope:           mean slope for each lap through firing field
%               lapR2:              mean R-Squared for each lap
%               lapCorrelation:     mean correlation for each lap (Kempter et al. 2012)
%               lapPhaseOffset:     mean phase offset foe each lap (Kempter et al. 2012)
%               circLinCorr:        Correlation between position and phase (Kempter et al. 2012)
%               pval:               p value for circ lin correlation (Kempter et al. 2012)
%               slopeCpU:           slope of regression line in degrees calculated by Kempter method (Kempter et al. 2012)
%               phaseOffset:        phase offset (Kempter et al. 2012)
%               RR:                 R-length (Kempter et al. 2012)

% ----------------------------------------------------------------------------------------
% Ryan E Harvey April 2017; 
% edited April 25th 2018; 
% edited Dec 2nd 2018: to accept multiple fields
%
% ADD FMA TO PATH
com=which('PHprecession');
com=strsplit(com,filesep);
basedir=[com{1},filesep,'Users',filesep,com{3},filesep,'GoogleDrive',filesep,'MatlabDir'];
addpath(genpath([basedir,filesep, 'buzcode', filesep, 'externalPackages', filesep, 'FMAToolbox']))

ThPrecess.phaselock.Rlength=NaN;
ThPrecess.phaselock.Pval=NaN;
ThPrecess.slope=NaN;
ThPrecess.RSquared=NaN;
ThPrecess.scatteredPH=NaN;
ThPrecess.lapSlope=NaN;
ThPrecess.lapR2=NaN;
ThPrecess.lapCorrelation=NaN;
ThPrecess.lapPhaseOffset=NaN;
ThPrecess.circLinCorr=NaN;
ThPrecess.pval=NaN;
ThPrecess.slopeCpU=NaN;
ThPrecess.phaseOffset=NaN;
ThPrecess.stats=NaN;
ThPrecess.data=NaN;



if size(spks_VEL,1)<10
    return
end

% NORM POSITION
position=[occ4Ph(:,1),rescale(occ4Ph(:,2),0,1),rescale(occ4Ph(:,3),0,1)];

% CHECK TO SEE IF AT LEAST 10 SPIKES ARE WITHIN THE FIELD BOUNDARIES
[ts,idx]=unique(position(:,1));
rescalexspk=interp1(ts,position(idx,2),spks_VEL(:,1));
if sum(rescalexspk>fieldbound(1) & rescalexspk<fieldbound(2))<10
    return
end

% RUN FMA PRECESSION CODE
[data,stats]=PhasePrecession(position(idx,:),spks_VEL(:,1),phase,'boundaries',fieldbound);
ThPrecess.stats=stats;
ThPrecess.data=data;


% PLOT
% figure;
% PlotPhasePrecession(data,stats)
% REMOVE FMA
rmpath(genpath([basedir,filesep, 'buzcode', filesep, 'externalPackages', filesep, 'FMAToolbox']))

spks_VEL_working = interp1(phase(:,1),phase(:,2),spks_VEL(:,1)','linear');

% COMPUTE PHASE LOCKING
ThPrecess.phaselock.Rlength=circ_r(spks_VEL_working');
[ThPrecess.phaselock.Pval,~]=circ_rtest(spks_VEL_working');

% PLACE OUTPUT IN STRUCTURE
ThPrecess.slope=stats.slope;
ThPrecess.RSquared=stats.r2;
ThPrecess.scatteredPH=[data.position.x,wrapTo360(rad2deg(data.position.phase))];
lapslope=stats.lap.slope;
lapslope(isnan(lapslope))=[];
lapslope(isinf(lapslope))=[];
ThPrecess.lapSlope=nanmean(lapslope);

lapr2=stats.lap.r2;lapr2(isnan(lapr2))=[];lapr2(isinf(lapr2))=[];
ThPrecess.lapR2=nanmean(lapr2);

circ_lin_corr=[];
phi0_deg=[];
for i=1:length(stats.lap.slope)
    if length(data.position.x(data.position.lap==i))>=2
        [circ_lin_corr(i,1),pval(i,1),slope_deg(i,1),phi0_deg(i,1)]=...
            kempter_lincirc(data.position.x(data.position.lap==i),...
            rad2deg(data.position.phase(data.position.lap==i)));
    else
        circ_lin_corr(i,1)=NaN;pval(i,1)=NaN;slope_deg(i,1)=NaN;phi0_deg(i,1)=NaN;RR(i,1)=NaN;
    end
end
ThPrecess.lapCorrelation=nanmean(circ_lin_corr);
ThPrecess.lapPhaseOffset=nanmean(phi0_deg);

[ThPrecess.circLinCorr,ThPrecess.pval,ThPrecess.slopeCpU,ThPrecess.phaseOffset]=kempter_lincirc(data.position.x,data.position.phase);
end
% OLD CODE
% figure
% ylim([0 360])
% for i=0:360
%     tempph=((data.position.phase*180/pi)+i);
%     tempph(tempph>360)=tempph(tempph>360)-360;
%     
% plot(data.position.x,tempph,'.k')
% [cors(i+1),p(i+1)]=corr(data.position.x,tempph);
% lsline
% pause(.01)
% end
% max(abs(cors))
% figure
% scatter(cors,p)
%
% Hippocampal spatio-temporal receptive fields were computed by dividing
% the number of spikes fired by a neuron in each spatio-phase bin 
% (spatial width 1 pixel=0.66 cm, phase width 1 degree) by the total time spent by the rat in that bin.
%     
%     pos_edge=linspace(min(occ4Ph(:,2)),max(occ4Ph(:,2)),120/.66);
%     ph_edge=0:360;
%     
%     ph_pos_spk=histcounts2(rad2deg(data.position.phase),data.position.x,ph_edge,pos_edge);
%     ph_pos_occ=histcounts2(interp1(phase(:,1),rad2deg(phase(:,2)),occ4Ph(:,1)),occ4Ph(:,2),ph_edge,pos_edge);
% 
%     spatio_phase=ph_pos_spk./ph_pos_occ;
% 
%     
%     
%     
%       filtWidth = [360 7]; filtSigma = 40;
%     imageFilter=fspecial('gaussian',filtWidth,filtSigma);
%     
%     N = nanconv([spatio_phase;spatio_phase;spatio_phase],imageFilter, 'nanout');
%     N=N(721:1440,:);
%     
% 
%     % The result was smoothed by convolution with a two-dimensional gaussian (spatial width 2.4 cm, temporal width 12 degrees).
%     filtWidth = [12,4]; filtSigma = 100;
%     imageFilter=fspecial('gaussian',filtWidth,filtSigma);
%     spatio_phase = nanconv(spatio_phase,imageFilter, 'nanout');


% PLOT FOR DEBUGGING
% figure;plot(position(:,2),'.k')
% hold on
% scatter(interp1(position(:,1),[1:length(position)],spks_VEL(:,1)),interp1(position(:,1),position(:,2),spks_VEL(:,1)),'r')
% plot([1;length(position)],[fieldbound(2),fieldbound(2)])
% plot([1;length(position)],[fieldbound(1),fieldbound(1)])

%spikes on path
% figure;plot(position(:,2),position(:,3),'.k');hold on;scatter(interp1(position(:,1),position(:,2),spks_VEL(:,1)),interp1(position(:,1),position(:,3),spks_VEL(:,1)),'r')

% figure;
% for i=1:length(stats.lap.slope  )
%     x=data.position.x(data.position.lap==i);
% 
%     subplot(1,2,1)
%     scatter(x,data.position.phase(data.position.lap==i),'k');hold on
% 
%     len = max(x)-min(x);
%     x = mean([min(x);max(x)]);
%     y = wrapTo2Pi(stats.lap.slope(i)*x+stats.lap.intercept(i));%*180/pi;
% 
%     X = x+[-len len]/2;
%     Y = y+[-len len]/2*stats.lap.slope(i);
%     plot(X,Y,'r');
%     hold off
% 
%     subplot(1,2,2)
%     plot(1:i,stats.lap.slope(1:i),'k');hold on
%     pause(.1)
% 
%     if length(data.position.x(data.position.lap==i))>=2
%         [circ_lin_corr(i,1),pval(i,1),slope_deg(i,1),phi0_deg(i,1),RR(i,1)]=...
%             cl_corr(data.position.x(data.position.lap==i),...
%             rad2deg(data.position.phase(data.position.lap==i)),...
%             deg2rad(min(stats.lap.slope)),deg2rad(max(stats.lap.slope)));
%     else
%        circ_lin_corr(i,1)=0;pval(i,1)=NaN;slope_deg(i,1)=0;phi0_deg(i,1)=NaN;RR(i,1)=0;
%     end
% end


% WRAP ANGLES
% WrappedPH=wrapTo360(rad2deg(spks_VEL_working)');

% combine distance and phase for use outside function
% scatteredPH=[normalizedD, WrappedPH];

% Circular-Linear Correlation

% figure;
% RSquared=[];
% corrs=[];
% for i=0:359
%     tempphase=WrappedPH+i;
%     tempphase(tempphase>360)=tempphase(tempphase>360)-360;
%     corrs=[corrs;corr2(normalizedD,tempphase)];
%
% mdl=fitlm(normalizedD,tempphase);
% RSquared=[RSquared;mdl.Rsquared.Ordinary];
%
% % plot(normalizedD,tempphase,'.k'); hold on
% % plot(normalizedD,mdl.Fitted,'r')
% % pause(.1)
% % hold off
% end

% [m,I]=max(abs(corrs));
% tempphase=WrappedPH+I;
% tempphase(tempphase>360)=tempphase(tempphase>360)-360;
% % Correlation=corr2(normalizedD,tempphase);
% [R,P] = corrcoef(normalizedD,tempphase);
% Correlation=R(1,2);
% sigCorr=P(1,2);
% mdl=fitlm(normalizedD,tempphase);
% RSquared=mdl.Rsquared.Ordinary;
% slope = normalizedD\tempphase;
% slope=mdl.Coefficients.Estimate(2);

% figure
% plot(normalizedD,tempphase,'.k'); hold on
% plot(normalizedD,mdl.Fitted)
% %
% figure;
% plot(normalizedD,tempphase,'.k'); hold on
% plot(normalizedD,tempphase+360,'.k'); hold on


% 	[beta,r2,p] = CircularRegression(normalizedD,spks_VEL_working);
% 	slope = beta(1);
% 	stats.intercept = beta(2) + stats.slope*dx; % Correction for circular/circular regression (see above)
% 	stats.r2 = r2;
% 	stats.p = p;

% figure;
% plot(data.position.x,data.position.phase ,'.k' );hold on
% plot(data.position.x,data.position.phase+2*pi ,'.k' );hold off


% save coefficient of determination (R-Squared) values for future use
% [RSquared,I]=max(RSquared);
% mdl=fitlm(scatteredPH(:,1),circshift(scatteredPH(:,2),I));
% Correlation=corr2(circshift(scatteredPH(:,2),I),scatteredPH(:,1));

% obtain slope
% slope = scatteredPH(:,1)\circshift(scatteredPH(:,2),I);
% SMOOTHED MAP
% nBinsx = round(track_length/5);
% Smatrix = hist3([WrappedPH,spks_VEL(:,2)],[15,nBinsx]);
% % Calc Depth of Modulation (Terrazas et al., 2005)
% % minNonZero=sum(Smatrix,2); minNonZero(minNonZero<1)=[];
% % DOM=(max(sum(Smatrix,2))-min(minNonZero))/max(sum(Smatrix,2));
% % remove nans or infs
% Smatrix(isnan(Smatrix)) = 0; Smatrix(isinf(Smatrix)) = 0;
% % smooth
% filtWidth = [5 3]; filtSigma = 1; imageFilter=fspecial('gaussian',filtWidth,filtSigma);
% smoothedPHmap = nanconv(Smatrix,imageFilter, 'nanout');
% % get mean firing rate
% % meanFR=mean2(smoothedPHmap);
% % conbine maps to get 0:720 degrees
% smoothedPHmap=[smoothedPHmap;smoothedPHmap];

% % PLACE OUTPUT IN STRUCTURE
% ThPrecess.Correlation=Correlation;
% ThPrecess.sigCorr=sigCorr;
% ThPrecess.slope=slope;
% ThPrecess.smoothedPHmap=smoothedPHmap;
% ThPrecess.mdl=mdl;
% ThPrecess.RSquared=RSquared;
% ThPrecess.scatteredPH=scatteredPH;

% [ ThetaStats ] = ThetaPower(EEG_DownSampledData);

