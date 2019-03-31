function [ mean_vector_length,peak_Firing_Rate,preferred_Direction,halfPeak,Directional_Range_HalfWidth,Direct_infoContent,BinsNbSpikes,BinsAngle3,BinsAngle] = HDCell(spks,angle,sampleRate)
%HDCell Summary of this function goes here
 % extract spike directions from smoothed and velocity filtered data

        % number of directional bins
        dBins = 60; 
        
        % Binning of directional data can also be done using histc function
        for j = 1:dBins;
            k = 0:dBins;
            % number of instances where the head of the rat was found in each of the 60 possible orientations
            ListOrientation = find(angle >= k(j)*6 & angle < (k(j)+1)*6);
            if length(ListOrientation) < 1; % if the number is less than 1 (i.e., 0), it is hard set to 1 to avoid division by 0
                nOrientation(j) = 1;
                nSpikesOrientation(j) = 0; % 0 is assigned to the number of spikes for this orientation
            else
                nSpikesOrientation(j) = sum(spks(ListOrientation));
                nOrientation(j) = length(ListOrientation);
            end
        end
        
        % transformed values from 1/60th of a sec to seconds
        nOrientation2 = nOrientation./sampleRate; % 30ms for the time in s
        
        % calculates the spikes/sec for each 6° bin
        BinsNbSpikes = nSpikesOrientation./nOrientation2;
        
        % trace the graph of the discharge rate / direction to the cells 1
        BinsAngle = 0.052358333:0.104716667:6.283;
        BinsAngle3 = 3:6:360;
        
        % calculate head direction cell properties
        mean_vector_length = circ_r(BinsAngle',BinsNbSpikes',circ_ang2rad(6)); % mean vector length based on binned firing rates
        peak_Firing_Rate = max(BinsNbSpikes); % peak firing rate
        pfdi = find(BinsNbSpikes(1,:) == peak_Firing_Rate);
        preferred_Direction = BinsAngle3(pfdi); % preferred firing direction
        halfPeak = peak_Firing_Rate/2;
        hpi = find(BinsNbSpikes(1,:) >= halfPeak);
        Directional_Range_HalfWidth_bins = BinsAngle3(hpi);
        Directional_Range_HalfWidth = max(Directional_Range_HalfWidth_bins) - min(Directional_Range_HalfWidth_bins);
        
        % calculate directional information content (from Taube & Muller 1998)
        probOrient = nOrientation2./sum(nOrientation2); % probability of occupancy
        overall_rate = (sum(nSpikesOrientation))/(sum(nOrientation2));
        reIC = BinsNbSpikes'./overall_rate;
        log_IC = log2(reIC);
        ij= find(isinf(log_IC)); % find -Inf's (log(0)) and replace with 0's (based on code from McN lab)
        log_IC(ij) = 0;
        ICi = probOrient'.*reIC.*log_IC;
        Direct_infoContent = sum(ICi);

end

