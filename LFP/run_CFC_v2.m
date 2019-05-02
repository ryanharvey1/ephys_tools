%Cross-Frequency coupling script for B.Clark lab data structure
%
%

%LB October 2018

function cfc=run_CFC_v2(data)
% run_CFC will compute max MI for spectral analysis and MI for two
% a priori questions about specific bands via tort 2010. Code adapted from
% CMBHome toolbox.
%
%
% Dependencies: getPower.m, thetaModGamma_v2.m, UnderThreshold Detect_v2.m
data=load(data);

%ADD PATHS FOR PROCESSING
com=which('run_CFC');
com=strsplit(com,filesep);

basedir=[com{1},filesep,'Users',filesep,com{3},filesep,'GoogleDrive',filesep,'MatlabDir'];
addpath([basedir,filesep,'BClarkToolbox',filesep,'Analysis'],...
    [basedir,filesep,'BClarkToolbox',filesep,'Analysis',filesep,'Utils'],...
    [basedir,filesep,'BClarkToolbox',filesep,'Analysis',filesep,'Visualize'],...
    [basedir,filesep,'BClarkToolbox',filesep,'Analysis',filesep,'LFP']);
addpath(genpath([basedir,filesep,'chronux_2_11']))
clear com

%SET PARAMETERS
% nchannels=1:size(data.lfp.signal,1);
nsessions=data.events;
mazeType=data.mazetypes;

%Calculate CFC per session
for ses=1
    
    %PARSE DATA FOR SESSION USING EVENTS
    session.mazetype=mazeType{1,ses};
    
    session.signal=data.lfp.signal(:,data.lfp.ts>=data.events(1,ses)...
        & data.lfp.ts<=data.events(2,ses));
    
    session.ts=data.lfp.ts(:,data.lfp.ts>=data.events(1,ses)...
        & data.lfp.ts<=data.events(2,ses));
    
    session.lfpsamplerate=data.lfp.lfpsamplerate;
    
    session.vel=data.frames(data.frames(:,1)>=data.events(1,ses)...
        & data.frames(:,1)<=data.events(2,ses),5);
    
    session.tsVel=data.frames(data.frames(:,1)>=data.events(1,ses)...
        & data.frames(:,1)<=data.events(2,ses),1);
    
    %FIND MOST ACTIVE CHANNEL FOR THETA
    %     disp('Finding most active channel for theta')
    [session.theta_ch.f_peak, session.theta_ch.active_lfp_ind] = ChooseActiveLFP_v2(session);
    disp(['Got it! Channel ',num2str(session.theta_ch.active_lfp_ind),' is the best'])
    
    % Get whole session MI (Spectral Analysis)
    
    session.signal=session.signal(session.theta_ch.active_lfp_ind,:);
    % RUN BROAD BAND CFC ON MULTIPLE FREQUENCY BANDS FOR EXPLORATORY
    % ANALYSIS OF PHASE-AMPLITUDE COUPLING
    session.CFC = thetaModGamma_v2(session,...
        'filtType','vel','filtParams',[5 90],'thetarange',...
        4:0.25:12,'gammarange',25:1:55,'nbins',36,'ifPlot',0);
    
end

%     %Get rid of redundant data
remove={'signal','ts','lfpsamplerate','vel','tsVel'};
session=rmfield(session,remove);

%Save session data
cfc{ses}=session;

clear session

end
%
% % save mat file to temp file
% if contains(data.session_path,'HPCatn')
%     save(['D:\Projects\HPCatn\CFCResults\',data.rat,'_',data.sessionID],...
%         'cfc','-v7.3')
% elseif contains(data.session_path,'PAE')
%     save(['D:\Projects\PAE_PlaceCell\CFCResults\',data.rat,'_',data.sessionID],...
%         'cfc','-v7.3')
% elseif contains(data.session_path,'ClarkP30')
%     save(['D:\ClarkP30_Recordings\CFCResults\',data.rat,'_',data.sessionID],...
%         'cfc','-v7.3')
% end

%######################## LOCAL FUNCTIONS ################################

%
% function [r_pp] = phase_phase(session,ch)
%
% low_filtered = eegfilt(session.signal(ch,:),session.lfpsamplerate,4,12);
% high_filtered = eegfilt(session.signal(ch,:),session.lfpsamplerate,12,16);
%
% phase1 = angle(hilbert(low_filtered));
% phase2 = angle(hilbert(high_filtered));
%
% r_pp = [];
% for sample = 1:300
%
%     upper_range = length(phase1) - 1*session.lfpsamplerate;
%     lower_range = 1;
%
%     I_rand = (upper_range-lower_range).*rand(1,1) + lower_range;
%     timewindow = round(I_rand):round(I_rand+1*session.lfpsamplerate);
%
%     r_sample = [];
%     counter = 0;
%     for mm = 1:20 % DataInput.par.nmcurve               = 1:20;
%         counter = counter + 1;
%         r_sample(:,counter) = phase_phase_core(phase1(timewindow),phase2(timewindow),mm);
%     end
%
%     r_pp(sample,:) = abs(r_sample);
% end
%
% end
%
% function [PhaseCount,binn] = phase_phase_histogram(session,ch)
%
% low_filtered = eegfilt(session.signal(ch,:),session.lfpsamplerate,6,12);
% high_filtered = eegfilt(session.signal(ch,:),session.lfpsamplerate,30,90);
%
% phase1 = angle(hilbert(low_filtered));
% phase2 = angle(hilbert(high_filtered));
%
% binwidth = pi/60;
% phasebin = -pi:binwidth:(pi-(binwidth));
% PhaseCounts = [];
% for j = 1:length(phasebin)
%     I=find(phase1 >= phasebin(j) & phase1 < phasebin(j)+binwidth);
%     [counts, ~] = hist(phase2(I), phasebin+binwidth/2);
%     PhaseCounts(j,:)=counts;
%
% end
%
% binn = phasebin+binwidth/2;
% end
%
% function plot_phase_histogram(binn,PhaseCounts )
% phasebin = binn;
% Phases = PhaseCounts ;
%
% figure
% myfilter = fspecial('gaussian',[10 10], 20);
% myfilteredimage = imfilter(Phases', myfilter, 'replicate');
% imagesc(phasebin,phasebin,myfilteredimage)
% axis xy square
% colorbar
% set(gcf,'color','w')
% xlabel('Slower Oscillation Phases')
% ylabel('Faster Oscillation Phases')
% title('Phase Counts')
%
% end

%% CODE GRAVEYARD


%     %DEFINE EPOCHS WHERE THETA DELAT RATIO IS ABOVE THRESHOLD (e.g. rat is
%     %moving)
%     disp('Find epochs with good theta delta ratio')
%     [session.epochs, session.P_r, session.epochTs]=ThetaEpochs_v2(session,active_lfp_ind);

%compute MI for theta-phase low-gamma amplitude coupling
%     [session.CFC.LGmodindex, session.CFC.LGmeanAmps, session.CFC.LGthetarange,...
%         session.CFC.LGgammarange,session.CFC.LGpartial] = ...
%         thetaModGamma_v2(session,active_lfp_ind,'filtType','vel','filtParams',[5 90],'thetarange',...
%         4:0.25:12,'gammarange',30:1:55,'ifPlot',1);
%
%     %compute MI for theta-phase high-gamma amplitude coupling
%         [session.CFC.HGmodindex, session.CFC.HGmeanAmps, session.CFC.HGthetarange,...
%         session.CFC.HGgammarange,session.CFC.HGpartial] = ...
%         thetaModGamma_v2(session,active_lfp_ind,'filtType','vel','filtParams',[5 90],'thetarange',...
%         4:0.25:12,'gammarange',60:1:90,'ifPlot',1);

%RUN CFC CODE ON SURROGATE DATA (VIA SHIFTING THETA ANGLES)
%compute MI surrogate for theta-phase low-gamma amplitude coupling
%     [session.surrogateCFC.LGmodindex, session.surrogateCFC.LGmeanAmps,...
%         session.surrogateCFC.LGthetarange, session.surrogateCFC.LGgammarange,...
%         session.surrogateCFC.LGpartial] = ...
%         thetaModGamma_v2(session,active_lfp_ind,'filtType','vel','filtParams',[5 90],'shuffles'...
%         ,200,'thetarange',4:0.25:12,'gammarange',30:1:55,'ifPlot',1);
%
%     %compute MI surrogate for theta-phase high-gamma amplitude coupling
%         [session.surrogateCFC.LGmodindex, session.surrogateCFC.LGmeanAmps,...
%         session.surrogateCFC.LGthetarange, session.surrogateCFC.LGgammarange,...
%         session.surrogateCFC.LGpartial] = ...
%         thetaModGamma_v2(session,active_lfp_ind,'filtType','vel','filtParams',[5 90],'shuffles'...
%         ,200,'thetarange',4:0.25:12,'gammarange',60:1:90,'ifPlot',1);
%
%         RUN CFC CODE ON SURROGATE DATA (VIA SHIFTING THETA ANGLES)
%       compute MI surrogate for theta-phase low-gamma amplitude coupling
% [session.surrogateCFC.modindex, session.surrogateCFC.meanAmps,...
%     session.surrogateCFC.thetarange, session.surrogateCFC.gammarange,...
%     session.surrogateCFC.partial] = ...
%     thetaModGamma_v2(session,session.theta_ch.active_lfp_ind,'filtType','vel','filtParams',[5 90],'shuffles'...
%     ,200,'thetarange',4:0.25:12,'gammarange',30:1:200,'nbins',36,'ifPlot',1);