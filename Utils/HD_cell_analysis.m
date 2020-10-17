classdef HD_cell_analysis
    % Class of functions to quantify features of head direction tuning
    
    % To do: 
    %   Measures to create: 
    %       - Anticipatory Time Interval
    %       - Tuning Quality (Preferred Direction(s),tuning width, kappa)
    %       - Ensemble Analysis 
    % Misc: 
    %       - Standardize input arguments (data,session,cell)
    
    properties
        Property1
    end
    
    methods(Static)
        
   % Distributive Ratio
    
        function DR = distributiveRatio(ratemap,frames,hdTuning,mazesize)
            
            thetabins = rad2deg((pi/30)/2:(pi/30):2*pi-(pi/30)/2);
            
            nBinsx = round(mazesize/3); nBinsy = round(mazesize/3);
            MinY = min(frames(:,3));
            MaxY = max(frames(:,3));
            MinX = min(frames(:,2));
            MaxX = max(frames(:,2));
            edges{1} = linspace(MinY, MaxY, nBinsy+1);
            edges{2} = linspace(MinX, MaxX, nBinsx+1);
            
            % Build
            for i = 1:length(thetabins)-1
                % find frames in specific theta bin
                frames_by_theta = frames(frames(:,4) > thetabins(i) & frames(:,4) < thetabins(i+1),:);
                % bin xy
                theta_map = hist3([frames_by_theta(:,3) frames_by_theta(:,2)],'Edges',edges);
                theta_map(end,:) = [];
                theta_map(:,end) = [];
                
                PR(i) = nansum(ratemap.*theta_map)/nansum(theta_map(:)); %Calculates predicted rate
            end
            
            for i = 1:length(hdTuning)
                tempDR(i) = log((1+hdTuning(i)))/(1+PR(i));
            end
            
            DR = nansum(abs(tempDR))/length(hdTuning);
            
        end
        
   % Directional information content
     
        function DIC  = computeDIC(histAng,hdTuning,OverallFR)
            %Computes the directional information content on unsmoothed ratemap for
            %directional data. Adapted from Langston et al 2010. by LB March 2018
            %
            %Input:
            %    histAng:   histogram of time spent in angular tuning
            %               curve.
            %    hdTuning:  head direction tuning curve.
            %    OverallFR: overall mean firing rate of the cell.
            %Output:
            %   DIC: directional information content in bits/spk
            %
            % by LB March 2018
            
            probOrient = histAng./sum(histAng);
            reIC = hdTuning./OverallFR;
            log_IC = log2(reIC); log_IC(isinf(log_IC)) = 0;
            ICi = probOrient.*reIC.*log_IC;
            DIC = sum(ICi)/OverallFR;
        end
        
   % Tuning Stability
    
        function stable_score = stability(data_video_spk,samplerate)
            % Computes the correlation of tuning curves generated after each
            % complete sampling of the horizontal azimuth.
            
            %Input:
            %     data_video_spk: n x 6 matrix containing timestamps, x
            %     coords, y coords, head direction (deg), instantaneous
            %     velocity (cm/s), and spike binary. From B.Clark lab data
            %     structure (obtained using createframes_w_spikebinary).
            %
            %Output:
            %   stability: mean correlation between all pairwise
            %   correlations of tuning curves
            %
            % By LB & RH 2018, updated by LB September 2019, Updated by RH
            % 2020
            
            %initlize loop parameters
            a = data_video_spk(data_video_spk(:,6) == 0,4);
            ts = data_video_spk(data_video_spk(:,6) == 0,1);
            frames = data_video_spk(data_video_spk(:,6) == 0,:);
          
            bins = 0:6:360;
            bin_check = zeros(1,length(bins));
            ii = 1;
            hd_tune = [];
            for i = 1:size(a,1)
                for b = 1:length(bins)
                    if a(i) > bins(b) && a(i) <= bins(b+1)
                        bin_check(b) = 1;
                        break
                    end
                end
                if sum(bin_check) == length(bin_check)-1
                    idx = frames(:,1) >= ts(ii) & frames(:,1) <= ts(i);
                    temp_a = a(idx);
                    
                    idx = data_video_spk(:,1) >= ts(ii) & data_video_spk(:,1) <= ts(i);
                    temp_spk_frames = data_video_spk(idx,:);
                    
                    [~,~,~,~,~,hdtuning] = tuningcurve(temp_a,...
                        temp_spk_frames(temp_spk_frames(:,6) == 1,4),samplerate);
                    
                 
                    hd_tune = [hd_tune;hdtuning];
                    
                    ii = i+1;
                    bin_check = zeros(1,length(bins));
                end
                if i == size(a,1) && isempty(hd_tune)
                    stable_score = NaN;
                    warning('insufficient azimuth sampling')
                    return
                end
            end
            
            %smooth over 36 degrees
            hd_tune = smoothdata([hd_tune hd_tune hd_tune],2,'gaussian',6);
            hd_tune = hd_tune(:,61:120);
            
            corrMat = corr(hd_tune');
            corrMat(logical(eye(length(corrMat)))) = NaN;
            stable_score = nanmean(corrMat(:));
        end
        
   % Tuning Stability across session
   
        function [within_Coeff,within,normWithin] = four_quarter_stability(data_video_spk,sampleRate,method)
            %four_quarter_stability computes the 4-quarter stability score for head direction signals based off
            %Boccara et al.
            
            %INPUT:
            %       -data_video_spk: timestamps interpolated with spike timestamps
            %       -sampleRate: video sample rate (in hz)
            %       -Angle:
            %       -Spike:
            
            %OUTPUT:
            %       -within_Coeff: 4-qrt stability score
            %       -structure of raw tuning curves for each qtr in ascending order
            %       (e.g. row 1= qtr 1).
            %       -normWithin: matrix of normalize tuning curves (row 1= qtr 1, row
            %       2=qtr 2 etc.)
            
            % Created by LBerkowitz March 2018, updated by LB July 2018, updated by LB fixed bug
            % August 2019
            
            time = [.25,.5,.75,1];
            
            start = min(data_video_spk(:,1));
            stop = max(data_video_spk(:,1));
            blocking_factor = stop-start;
            
            block = blocking_factor*time;
            
            clear blocking_factor start stop time
            
            block = [data_video_spk(1,1) data_video_spk(1,1)+block];
            
            for i = 1:length(block)-1
                tempSpk = data_video_spk(data_video_spk(:,1) > block(i) & data_video_spk(:,1) < block(i+1),:);
                
                [~,~,~,~,~,hdTuning] = tuningcurve(tempSpk(tempSpk(:,6)==0,4),tempSpk(tempSpk(:,6)==1,4),sampleRate);
                
                within.hdTuning{i,1} = hdTuning;
                
            end
            clear block
            
            temp = [within.hdTuning{1,1}; within.hdTuning{2,1}; within.hdTuning{3,1}; within.hdTuning{4,1}];
            
            normWithin = (temp - min(temp(:))) / (max(temp(:)) - min(temp(:)));
            
            if strcmp(method,'cor')
                
                %calculate correlation for all pairs
                first=corr2(within.hdTuning{1,1},within.hdTuning{2,1});
                second=corr2(within.hdTuning{1,1},within.hdTuning{3,1});
                third=corr2(within.hdTuning{1,1},within.hdTuning{4,1});
                fourth=corr2(within.hdTuning{2,1},within.hdTuning{3,1});
                fifth=corr2(within.hdTuning{2,1},within.hdTuning{4,1});
                sixth=corr2(within.hdTuning{3,1},within.hdTuning{4,1});
                
                %warning if data is nan
                if sum(isnan([first,second,third,fourth,fifth,sixth]))>1
                    warning('NaNs within the data are preventing an accurate calculation of stability. Check your data')
                    
                end
                
                within_Coeff=nanmean([first,second,third,fourth,fifth,sixth]);
                
            elseif strcmp(method, 'std')
      
                within_Coeff = nanmean(nanstd(normWithin));
                
            end
            
%             normTemp=[];
%             for i=1:4
%                 tempTune=within.hdTuning{i,:};
%                 normTemp=[normTemp; rescale(tempTune,0,1)];
%             end
%             normWithin=normTemp;
        end
   
   % Anticipatory time interval (NOT COMPLETE)
   
        function anticipatory_time_interval = computeATI(data,session,cell)
            
            
            [data_video_spk,data_video_nospk]=createframes_w_spikebinary(data,session,cell);
            
            %Compute angular head velocity
            anglevel=insta_angvel(data_video_nospk(:,4),30);
            
            %Interpolate head direction to insure HD is in phase with
            %velocity
            new_time = data_video_nospk(:,1) + ((1/data.samplerate)/2); %set time inbetween data points
            new_time(new_time>max(data_video_nospk(:,1)))=[];
            
            in_phase_hd = interp1(data_video_nospk(:,1),data_video_nospk(:,4),new_time,'linear');
            in_phase_ang_vel = interp1(new_time,anglevel,new_time,'linear');
            
            
            spk_ts = data_video_spk(data_video_spk(:,6)==1,1);
            
            in_phase_vel_spk = interp1(new_time,anglevel,...
                spk_ts,'linear');
            in_phase_hd_spk = interp1(new_time,in_phase_hd,...
                spk_ts,'linear');
            
            % Identify clockwise (v > 90deg/s), counterclockwise(v > -90deg/s), and still frames |v| < 30deg/s and
            % make tuning curves.
            
            [~,~,~,~,prefdirec,hdTuning_c] =...
                tuningcurve(in_phase_hd(anglevel > 90),in_phase_hd_spk(in_phase_vel_spk > 90),30)
            
            [~,~,~,~,prefdirec,hdTuning_cc] =...
                tuningcurve(in_phase_hd(anglevel > 90),in_phase_hd_spk(in_phase_vel_spk < -90),30)
            
            [~,~,~,~,prefdirec,hdTuning_s] =...
                tuningcurve(in_phase_hd(anglevel > 90),in_phase_hd_spk(abs(in_phase_vel_spk) < 30),30)
            
            [~,~,~,~,prefdirec,hdTuning] =...
                tuningcurve(data_video_nospk(:,4),data_video_spk(data_video_spk(:,6)==1,4),30)
            
            % ADD CROSS_VALIDATION FOR ABOVE TO VERIFY MEAN DIRECTION IS
            % NOT OBTAINED BY CHANCE
            
            figure;
            plot(hdTuning_s,'Color',[.7 .7 .7],'LineWidth',2)
            hold on;
            plot(hdTuning_c,'b','LineWidth',2)
            plot(hdTuning_cc,'r','LineWidth',2)
            plot(hdTuning,'k','lineWidth',3)
            legend({'Still','clockwise','counter clockwise','Overall'})
            ylabel('Firing Rate')
            xlabel('Binned Direction')
            
            %Identify subcategories of fast and slow
            %Clkfast (v >= 270 deg/s)
            %Clkslow (30 deg/s <= v < 270 deg/s)
            %CCfast (v >= -270 deg/s)
            %CCslow (-30 deg/s <= v < -270 deg/s)
            %Still (-30 deg/s < v < 30 deg/s)
            
            %Compute angular difference between tuning curves of each
            %category
            
            %time-shift analysis for future (add 1/samplerate to ts) and past
            %(subtract 1/samplerate to ts) and current

        end
        
   
        
    end
end

