classdef HD_cell_analysis
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods(Static)
        
        %Distributive Ratio
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
        
        function stability=stability(data_video_spk,samplerate)
            %Computes the correlation of tuning curves generated after each
            %complete sampling of the horizontal azimuth.
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
            % By LB & RH 2018, updated by LB September 2019
            
            %initialize bins to determine sampling
            angBins = rad2deg((pi/30)/2:(pi/30):2*pi-(pi/30)/2);
            
            %initlize loop parameters
            i = 1;
            ii = 2;
            chunk = [];
            hdTuning = [];
            heading = data_video_spk(:,4);
            
            figure;
            plot(data_video_spk(:,2),data_video_spk(:,3),'.k')
            
            while ii < length(heading)
                histAng = histcounts(heading(i:ii),angBins);
                if sum(histAng>0)==length(angBins)-1
                    temp = data_video_spk(i:ii,:);
                    spk_a = temp(temp(:,6)==1,4);
                    a = temp(temp(:,6)==0,4);
                    [~,~,~,~,~,hdtuning] = tuningcurve(a,spk_a,samplerate);
                    hdTuning = [hdTuning;hdtuning];
                    chunk = [chunk;i,ii];
                    hold on
                    plot(data_video_spk(i:ii,2),data_video_spk(i:ii,3))
                    i = ii+1;
                    ii = ii+1;
                end
                ii = ii+1;
            end
            
            store = [];
            corrMat = corr(hdTuning');
            for k = 1:size(corrMat,2)-1
                data = diag(corrMat,k);
                store = [store;data];
            end
            
            stability = nanmean(store);
            
        end
        
        function [within_Coeff,within,normWithin] = four_quarter_stability(data_video_spk,sampleRate)
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
            
            normTemp=[];
            for i=1:4
                tempTune=within.hdTuning{i,:};
                normTemp=[normTemp; rescale(tempTune,0,1)];
            end
            normWithin=normTemp;
            within_Coeff=nanmean([first,second,third,fourth,fifth,sixth]);
        end
        
        function anticipatory_time_interval=computeATI(data,session,cell)
            
            %Get frames 
            event = data.events(session);
            frames = data.frames; 
            
            
            %Compute angular head velocity
            
            %Interpolate head direction to insure HD is in phase with
            
            %Identify clockwise (v > 90deg/s), counterclockwise(v > -90deg/s), and still frames |v| < 30deg/s.  
            
            %compute mean direction for each group 
            
            %Identify subcategories of fast and slow
                %Clkfast (v >= 270 deg/s)
                %Clkslow (30 deg/s <= v < 270 deg/s)
                %CCfast (v >= -270 deg/s)
                %CCslow (-30 deg/s <= v < -270 deg/s)
                %Still (-30 deg/s < v < 30 deg/s)
                
            %Compute angular difference between tuning curves of each
            %category 
            
            %time-shift analysis for future (add 1/samplerate to ts) and past
            %(subtract 1/samplerate to ts) and current (
            
            
        end
        
    end
end

