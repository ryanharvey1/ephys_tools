% batch_SpikesDirection creates summary of directional analyses
%
% Input: 
%       - raw data output from Labview
%
% Output:
%       - mean_vector_length = mean vector length for angle data
%       - preferred_Direction
%       - peak firing rate
%       - FirstSec_Corr = correlation of tuning curve between first and second half of session
%       - Direct_infoContent
%       - overall_rate
%       - FourQuart_MeanCorr
%       - Directional_Range_HalfWidth
% 
% this script requires Chronux/Circular stats toolbox and FindFiles function
% 
% Created by Ben C Jan 2014

% identify path to data files
path = '/Users/bjclark/Desktop/Dropbox/EC_analysis/__RawData_&_Summary/__RateMapData/__lesion_data/lg_mec';
ReadData = FindFiles('bc254_s10_TT3_u1.txt', 'StartingDirectory', path);

for i = 1:length(ReadData);
    
    data = importdata(ReadData{i}); % load xy data
    
        % video frame rate in Hz
    sampleRate = 60;

    % maze parameters for non-detects
    xmin = 0;
    xmax = 255;
    ymin = 0;
    ymax = 255;

    % min and max LED distance parameters
    minLED = 4;
    maxLED = 39;

    % number of directional bins
    dBins = 60;

    % extract coords, spikes, angle, and direction from ReadData output
    rx = data.data(:,2); % red x-coords
    ry = data.data(:,3); % red y-coords
    gx = data.data(:,4); % green x-coords
    gy = data.data(:,5); % green y-coords
    spks = data.data(:,6); % spikes
    angle = data.data(:,10); % angle in radians
    distance = data.data(:,11); % distance in pixels between LEDs
    datai = [rx,ry,gx,gy,spks,angle,distance]; % create array with all variables

    % find red LED non-detects
    dataFiltx = find(datai(:,1) > xmin & datai(:,1) < xmax & datai(:,2) < ymax & datai(:,1) > ymin);
    rFILT = datai(dataFiltx,:);

    % find green LED non-detects
    dataFiltxy = find(rFILT(:,3) > xmin & rFILT(:,3) < xmax & rFILT(:,4) < ymax & rFILT(:,4) > ymin);
    rgFILT = rFILT(dataFiltxy,:);

    % find Min and Max LED distance
    dataFiltxLED = find(rgFILT(:,7) > minLED & rgFILT(:,7) < maxLED);
    rgmmFILT = rgFILT(dataFiltxLED,:);

    % extract x y coords and smooth
    rxf = rgmmFILT(:,1);
    ryf = rgmmFILT(:,2);
    rxs = runline(rxf,5,1); % smooth with 10pt window and 1pt step (based on Valente et al 2007 PloS One)
    rys = runline(ryf,5,1);

    % extract spike locations
    spksf = rgmmFILT(:,5);
    spksi = find(spksf == 1); 
    spks_xs = rxs(spksi);
    spks_ys = rys(spksi);
    spks_xf = rxf(spksi);
    spks_yf = ryf(spksi);

    % extract angle and convert to degrees
    anglef = rgmmFILT(:,6);
    degHeading = circ_rad2ang(anglef);
    dSpks = degHeading(spksi);
    aSpks = anglef(spksi);
    timestamps = 1:length(degHeading);
    timestamps = timestamps';
    tSpks = timestamps(spksi);

    % Binning of directional data can also be done using histc function
    for j = 1:dBins;
        % number of instances where the head of the rat was found in each of the 60 possible orientations
        ListOrientation = find(anglef >= ((((2*pi)/dBins)/2)+(j-1)*2*pi/dBins) & anglef < ((((2*pi/dBins)/2)+(j)*2*pi/dBins)));
            if length(ListOrientation) < 1; % if the number is less than 1 (i.e., 0), it is hard set to 1 to avoid division by 0
                nOrientation(j) = 1; 
                nSpikesOrientation(j) = 0; % 0 is assined to the number of spikes for this orientation
            else
                nSpikesOrientation(j) = sum(spksf(ListOrientation));
                nOrientation(j) = length(ListOrientation);
            end
    end

    % transformed values from 1/60th of a sec to seconds
    nOrientation2 = nOrientation./sampleRate; % 60ms for the time in s

    % calculates the spikes/sec for each 6° bin
    BinsNbSpikes = nSpikesOrientation./nOrientation2;

    % trace the graph of the discharge rate / direction to the cells 1
    BinsAngle = (0+((2*pi/(dBins)/2)):2*pi/(dBins):(2*pi)-((2*pi/(dBins)/2)));
    BinsAngle3 = BinsAngle*dBins;

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
    ij = find(isinf(log_IC)); % find -Inf's (log(0)) and replace with 0's (based on code from McN lab)
    log_IC(ij) = 0;
    ICi = probOrient'.*reIC.*log_IC;
    Direct_infoContent = sum(ICi);

    % xcorr of tuning function between first and second half of session
    sessLength_samples = numel(timestamps);
    First_End = round(sessLength_samples/2);
    Sec_Start = First_End + 1;
    Sec_End = sessLength_samples;
    First_spksf = spksf(1:First_End);
    Sec_spksf = spksf(Sec_Start:Sec_End);
    First_anglef = anglef(1:First_End);
    Sec_anglef = anglef(Sec_Start:Sec_End);

    for k = 1:dBins;
        ListOrientation_First = find(First_anglef >= ((((2*pi)/dBins)/2)+(k-1)*2*pi/dBins) & First_anglef < ((((2*pi/dBins)/2)+(k)*2*pi/dBins)));
            if length(ListOrientation_First) < 1;
                nOrientation_First(k) = 1; 
                nSpikesOrientation_First(k) = 0;
            else
                nSpikesOrientation_First(k) = sum(First_spksf(ListOrientation_First));
                nOrientation_First(k) = length(ListOrientation_First); 
            end
    end

    nOrientation_First2 = nOrientation_First./sampleRate;
    BinsNbSpikes_First = nSpikesOrientation_First./nOrientation_First2;

    for l = 1:dBins;
        ListOrientation_Sec = find(Sec_anglef >= ((((2*pi)/dBins)/2)+(l-1)*2*pi/dBins) & Sec_anglef < ((((2*pi/dBins)/2)+(l)*2*pi/dBins)));
            if length(ListOrientation_Sec) < 1;
                nOrientation_Sec(l) = 1; 
                nSpikesOrientation_Sec(l) = 0;
            else
                nSpikesOrientation_Sec(l) = sum(Sec_spksf(ListOrientation_Sec));
                nOrientation_Sec(l) = length(ListOrientation_Sec); 
            end
    end

    nOrientation_Sec2 = nOrientation_Sec./sampleRate;
    BinsNbSpikes_Sec = nSpikesOrientation_Sec./nOrientation_Sec2;

    FirstSec_Corr = corrcoef(BinsNbSpikes_First',BinsNbSpikes_Sec'); % correlation between first half and second half of session
    FirstSec_Corr = FirstSec_Corr(1,2);

    % xcorr of tuning function between four quarters of session
    F_End = round(sessLength_samples/4);
    S_Start = F_End + 1;
    S_End = F_End*2;
    Third_Start = S_End + 1;
    Third_End = F_End*3;
    Fourth_Start = Third_End + 1;
    Fourth_End = sessLength_samples;

    F_spksf = spksf(1:F_End);
    S_spksf = spksf(S_Start:S_End);
    Third_spksf = spksf(Third_Start:Third_End);
    Fourth_spksf = spksf(Fourth_Start:Fourth_End);

    F_anglef = anglef(1:F_End);
    S_anglef = anglef(S_Start:S_End);
    Third_anglef = anglef(Third_Start:Third_End);
    Fourth_anglef = anglef(Fourth_Start:Fourth_End);

    % first quarter
    for m = 1:dBins;
        ListOrientation_F = find(F_anglef >= ((((2*pi)/dBins)/2)+(m-1)*2*pi/dBins) & F_anglef < ((((2*pi/dBins)/2)+(m)*2*pi/dBins)));
            if length(ListOrientation_F) < 1;
                nOrientation_F(m) = 1; 
                nSpikesOrientation_F(m) = 0;
            else
                nSpikesOrientation_F(m) = sum(F_spksf(ListOrientation_F));
                nOrientation_F(m) = length(ListOrientation_F); 
            end
    end

    nOrientation_F2 = nOrientation_F./sampleRate;
    BinsNbSpikes_F = nSpikesOrientation_F./nOrientation_F2;

    %  second quarter
    for n = 1:dBins;
        ListOrientation_S = find(S_anglef >= ((((2*pi)/dBins)/2)+(n-1)*2*pi/dBins) & S_anglef < ((((2*pi/dBins)/2)+(n)*2*pi/dBins)));
            if length(ListOrientation_S) < 1;
                nOrientation_S(n) = 1; 
                nSpikesOrientation_S(n) = 0;
            else
                nSpikesOrientation_S(n) = sum(S_spksf(ListOrientation_S));
                nOrientation_S(n) = length(ListOrientation_S); 
            end
    end

    nOrientation_S2 = nOrientation_S./sampleRate;
    BinsNbSpikes_S = nSpikesOrientation_S./nOrientation_S2;

    %  third quarter
    for p = 1:dBins;
        ListOrientation_T = find(Third_anglef >= ((((2*pi)/dBins)/2)+(p-1)*2*pi/dBins) & Third_anglef < ((((2*pi/dBins)/2)+(p)*2*pi/dBins)));
            if length(ListOrientation_T) < 1;
                nOrientation_T(p) = 1; 
                nSpikesOrientation_T(p) = 0;
            else
                nSpikesOrientation_T(p) = sum(Third_spksf(ListOrientation_T));
                nOrientation_T(p) = length(ListOrientation_T); 
            end
    end

    nOrientation_T2 = nOrientation_T./sampleRate;
    BinsNbSpikes_T = nSpikesOrientation_T./nOrientation_T2;

    %  fourth quarter
    for q = 1:dBins;
        ListOrientation_Fourth = find(Fourth_anglef >= ((((2*pi)/dBins)/2)+(q-1)*2*pi/dBins) & Fourth_anglef < ((((2*pi/dBins)/2)+(q)*2*pi/dBins)));
            if length(ListOrientation_Fourth) < 1;
                nOrientation_Fourth(q) = 1; 
                nSpikesOrientation_Fourth(q) = 0;
            else
                nSpikesOrientation_Fourth(q) = sum(Fourth_spksf(ListOrientation_Fourth));
                nOrientation_Fourth(q) = length(ListOrientation_Fourth); 
            end
    end

    nOrientation_Fourth2 = nOrientation_Fourth./sampleRate;
    BinsNbSpikes_Fourth = nSpikesOrientation_Fourth./nOrientation_Fourth2;

    FS_corr = corrcoef(BinsNbSpikes_F',BinsNbSpikes_S'); 
    FT_corr = corrcoef(BinsNbSpikes_F',BinsNbSpikes_T'); 
    FFourth_corr = corrcoef(BinsNbSpikes_F',BinsNbSpikes_Fourth'); 
    SThird_corr = corrcoef(BinsNbSpikes_S',BinsNbSpikes_T'); 
    SFourth_corr = corrcoef(BinsNbSpikes_S',BinsNbSpikes_Fourth');
    TFourth_corr = corrcoef(BinsNbSpikes_T',BinsNbSpikes_Fourth');

    FourQuart_MeanCorr = (FS_corr(1,2) + FT_corr(1,2) + FFourth_corr(1,2) + SThird_corr(1,2) + SFourth_corr(1,2) + TFourth_corr(1,2))/6;
    
    keep('path', 'ReadData', 'sampleRate', 'xmin', 'xmax', 'ymin', 'ymax', 'minLED', 'maxLED', 'dBins', 'i','mean_vector_length', 'peak_Firing_Rate', 'preferred_Direction', 'overall_rate', 'Direct_infoContent', 'FirstSec_Corr', 'FourQuart_MeanCorr', 'Directional_Range_HalfWidth', 'Directional_Range_HalfWidth_bins');
    [filepath, filename] = fileparts(ReadData{i});
    save([filepath filesep filename '_HD_Properties.mat']);
end