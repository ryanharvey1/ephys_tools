%Cross-Frequency coupling script for B.Clark lab data structure

%LB October 2018, updated January 2019

function cfc=run_CFC(data)
% run_CFC will compute max MI for spectral analysis and MI for two
% a priori questions about specific bands via tort 2010. Code adapted from
% CMBHome toolbox.
%
%
% Dependencies: getPower.m, thetaModGamma_v2.m, UnderThreshold Detect_v2.m

data=load(data);

%SET PARAMETERS
nsessions=data.events;
mazeType=data.mazetypes;

%Calculate CFC per session
for ses=1:size(nsessions,2)
    
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
    [session.theta_ch.f_peak, session.theta_ch.active_lfp_ind] = ChooseActiveLFP_v2(session);
    
    %compute MI for only two frequency ranges for surrogate data via tort 2010 (code adapted
    %from CMBHOME thetaModGamma).
    session.signal=session.signal(session.theta_ch.active_lfp_ind,:);

    session.LGtortCFC = ...
        computeMI_hilbert(session,'filtType','vel','filtParams',[5 90],'shuffles'...
        ,200,'thetarange',6:0.25:10,'gammarange',25:1:55,'nbins',36,'ifPlot',0);
    
    session.HGtortCFC = ...
        computeMI_hilbert(session,'filtType','vel','filtParams',[5 90],'shuffles'...
        ,200,'thetarange',6:0.25:10,'gammarange',65:1:200,'nbins',36,'ifPlot',0);
    
%     
%     %Get rid of redundant data
    remove={'signal','ts','lfpsamplerate','vel','tsVel'};
    session=rmfield(session,remove);
    
    %Save session data
    cfc{ses}=session;
    
    clear session
    
end

% save mat file to temp file
if contains(data.session_path,'HPCatn')
    save(['D:\Projects\HPCatn\CFCResults\',data.rat,'_',data.sessionID],...
        'cfc','-v7.3')
elseif contains(data.session_path,'PAE')
    save(['D:\Projects\PAE_PlaceCell\CFCResults\',data.rat,'_',data.sessionID],...
        'cfc','-v7.3')
elseif contains(data.session_path,'ClarkP30')
    save(['D:\ClarkP30_Recordings\CFCResults\',data.rat,'_',data.sessionID],...
        'cfc','-v7.3')
end
end


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

%     if contains(session.mazetype,["linear_track","CircularTrack","LinearTrack"])
%         
%         %ADD LAP THRESHOLD CRITERA HERE
%         if (size(data.linear_track.lapinfo.laps,2)/2)>=10
%             
%             tempsession.lfpsamplerate=session.lfpsamplerate;
%             tempsession.xlimits=[min(data.linear_track.right{1, 1}.dataspks(:,2)) max(data.linear_track.right{1, 1}.dataspks(:,2))];
%             
%             %get right signal power
%             for lap=1:size(data.linear_track.right{1, 1}.laps,1)
%                 
%                 tempsession.signal=session.signal(session.theta_ch.active_lfp_ind,...
%                     session.ts > data.linear_track.right{1, 1}.laps{lap}(1,1)...
%                     & session.ts < data.linear_track.right{1, 1}.laps{lap}(end,1));
%                 
%                 tempsession.ts=session.ts(session.ts > data.linear_track.right{1, 1}.laps{lap}(1,1)...
%                     & session.ts < data.linear_track.right{1, 1}.laps{lap}(end,1));
%                 
%                 tempsession.tsVel=data.linear_track.right{1, 1}.laps{lap}(data.linear_track.right{1, 1}.laps{lap}(:,6)==0,1);
%                 tempsession.vel=data.linear_track.right{1, 1}.laps{lap}(data.linear_track.right{1, 1}.laps{lap}(:,6)==0,5);
%                 tempsession.xcoord=data.linear_track.right{1, 1}.laps{lap}(data.linear_track.right{1, 1}.laps{lap}(:,6)==0,2);
%                 
%                 session.power.rightLap.LG{lap} = getPower(tempsession,...
%                     'filtType','vel','filtParams',[5 90],'thetarange',...
%                     4:0.25:12,'gammarange',25:1:55);
%                 session.power.rightLap.HG{lap} = getPower(tempsession,...
%                     'filtType','vel','filtParams',[5 90],'thetarange',...
%                     4:0.25:12,'gammarange',65:1:200);
%                 
%             end
%             
%             clear lap tempsession
%             
%             %Get left signal power
%             tempsession.lfpsamplerate=session.lfpsamplerate;
%             tempsession.xlimits=[min(data.linear_track.left{1, 1}.dataspks(:,2)) max(data.linear_track.left{1, 1}.dataspks(:,2))];
%             
%             for lap=1:size(data.linear_track.left{1, 1}.laps,1)
%                 tempsession.signal=session.signal(session.theta_ch.active_lfp_ind,session.ts > data.linear_track.left{1, 1}.laps{lap}(1,1)...
%                     & session.ts < data.linear_track.left{1, 1}.laps{lap}(end,1));
%                 
%                 tempsession.ts=session.ts(session.ts > data.linear_track.left{1, 1}.laps{lap}(1,1)...
%                     & session.ts < data.linear_track.left{1, 1}.laps{lap}(end,1));
%                 
%                 tempsession.tsVel=data.linear_track.left{1, 1}.laps{lap}(data.linear_track.left{1, 1}.laps{lap}(:,6)==0,1);
%                 tempsession.vel=data.linear_track.left{1, 1}.laps{lap}(data.linear_track.left{1, 1}.laps{lap}(:,6)==0,5);
%                 tempsession.xcoord=data.linear_track.left{1, 1}.laps{lap}(data.linear_track.left{1, 1}.laps{lap}(:,6)==0,2);
%                 
%                 session.power.leftLap.LG{lap} = getPower(tempsession,...
%                     'filtType','vel','filtParams',[5 90],'thetarange',...
%                     4:0.25:12, 'gammarange',25:1:55);
%                 session.power.leftLap.HG{lap} = getPower(tempsession,...
%                     'filtType','vel','filtParams',[5 90],'thetarange',...
%                     4:0.25:12, 'gammarange',65:1:200);
%             end
%             
%             clear lap tempsession
%         end
%         
%         
%         % Get whole session MI (Spectral Analysis)
%         
% %         session.signal=session.signal(session.theta_ch.active_lfp_ind,:);
% %         % RUN BROAD BAND CFC ON MULTIPLE FREQUENCY BANDS FOR EXPLORATORY
% %         % ANALYSIS OF PHASE-AMPLITUDE COUPLING
% %         session.CFC = thetaModGamma_v2(session,...
% %             'filtType','vel','filtParams',[5 90],'thetarange',...
% %             4:0.25:12,'gammarange',25:1:120,'nbins',36,'ifPlot',0);  
%     else
%         
% %         %Just get whole session MI
% %         
% %         session.signal=session.signal(session.theta_ch.active_lfp_ind,:);
% %         
% %         % RUN BROAD BAND CFC ON MULTIPLE FREQUENCY BANDS FOR EXPLORATORY
% %         % ANALYSIS OF PHASE-AMPLITUDE COUPLING
% % % 
% %         session.CFC = thetaModGamma_v2(session,...
% %             'filtType','vel','filtParams',[5 90],'thetarange',...
% %             4:0.25:12,'gammarange',25:1:120,'nbins',36,'ifPlot',0);
% %
%     end


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