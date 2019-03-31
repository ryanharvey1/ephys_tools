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
path = 'E:\Experiments\EC\HDC';
ReadData = FindFiles('*.read', 'StartingDirectory', path);

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
    dataFiltx = find(datai(:,1) ~= xmin & datai(:,2) ~= xmax);
    rFILT = datai(dataFiltx,:);
    dataFiltxErr = find(rFILT(:,1) ~= xmax & rFILT(:,2) ~= xmax);
    rFILTerr = rFILT(dataFiltxErr,:);

    % find green LED non-detects
    dataFiltxy = find(rFILTerr(:,3) ~= ymin & rFILTerr(:,4) ~= ymax);
    rgFILT = rFILTerr(dataFiltxy,:);
    dataFiltxyErr = find(rgFILT(:,3) ~= ymax & rgFILT(:,4) ~= ymax);
    rgFILTerr = rgFILT(dataFiltxyErr,:);
   
    % find Min and Max LED distance
    dataFiltxLED = find(rgFILTerr(:,7) > minLED & rgFILTerr(:,7) < maxLED);
    rgmmFILT = rgFILTerr(dataFiltxLED,:);

    % extract spike locations and extract angle and convert to degrees
    spksf = rgmmFILT(:,5);
    anglef = rgmmFILT(:,6);
    degHeading = circ_rad2ang(anglef);
        
    % Binning of directional data can also be done using histc function
    for j = 1:dBins;
        k = 0:dBins;
        % number of instances where the head of the rat was found in each of the 60 possible orientations
        ListOrientation = find(degHeading >= k(j)*6 & degHeading < (k(j)+1)*6);
            if length(ListOrientation) < 1; % if the number is less than 1 (i.e., 0), it is hard set to 1 to avoid division by 0
                nOrientation(j) = 1; 
                nSpikesOrientation(j) = 0; % 0 is assigned to the number of spikes for this orientation
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
    BinsAngle = 0.052358333:0.104716667:6.283;
    BinsAngle3 = 3:6:360;

    % calculate head direction cell properties
    mean_vector_length = circ_r(BinsAngle',BinsNbSpikes',circ_ang2rad(6)); % mean vector length based on binned firing rates
    peak_Firing_Rate = max(BinsNbSpikes); % peak firing rate
    pfdi = find(BinsNbSpikes(1,:) == peak_Firing_Rate); 
    preferred_Direction = BinsAngle3(pfdi); % preferred firing direction
    halfPeak = peak_Firing_Rate/2;
    hpi = find(BinsNbSpikes(1,:) >= halfPeak); 
    Directional_Range_HalfWidth = length(hpi)*6;

    % calculate directional information content (from Taube & Muller 1998)
    probOrient = nOrientation2./sum(nOrientation2); % probability of occupancy
    overall_rate = (sum(nSpikesOrientation))/(sum(nOrientation2));
    reIC = BinsNbSpikes'./overall_rate;
    log_IC = log2(reIC);
    ij = find(isinf(log_IC)); % find -Inf's (log(0)) and replace with 0's (based on code from McN lab)
    log_IC(ij) = 0;
    ICi = probOrient'.*reIC.*log_IC;
    Direct_infoContent = sum(ICi);

    % xcorr of tuning function between four quarters of session
    sessLength_samples = numel(degHeading);
    First_End = round(sessLength_samples/4);
    Second_Start = First_End + 1;
    Second_End = First_End*2;
    Third_Start = Second_End + 1;
    Third_End = First_End*3;
    Fourth_Start = Third_End + 1;
    Fourth_End = sessLength_samples;

    First_spksf = spksf(1:First_End);
    Second_spksf = spksf(Second_Start:Second_End);
    Third_spksf = spksf(Third_Start:Third_End);
    Fourth_spksf = spksf(Fourth_Start:Fourth_End);

    First_degf = degHeading(1:First_End);
    Second_degf = degHeading(Second_Start:Second_End);
    Third_degf = degHeading(Third_Start:Third_End);
    Fourth_degf = degHeading(Fourth_Start:Fourth_End);
    
    % first quarter
    for j = 1:dBins;
        k = 0:dBins;
        ListOrientation_First = find(First_degf >= k(j)*6 & First_degf < (k(j)+1)*6);
            if length(ListOrientation_First) < 1;
                nOrientation_First(j) = 1; 
                nSpikesOrientation_First(j) = 0;
            else
                nSpikesOrientation_First(j) = sum(First_spksf(ListOrientation_First));
                nOrientation_First(j) = length(ListOrientation_First); 
            end
    end

    nOrientation_First2 = nOrientation_First./sampleRate;
    BinsNbSpikes_First = nSpikesOrientation_First./nOrientation_First2;

    %  second quarter
    for j = 1:dBins;
        k = 0:dBins;
        ListOrientation_Second = find(Second_degf >= k(j)*6 & Second_degf < (k(j)+1)*6);
            if length(ListOrientation_Second) < 1;
                nOrientation_Second(j) = 1; 
                nSpikesOrientation_Second(j) = 0;
            else
                nSpikesOrientation_Second(j) = sum(Second_spksf(ListOrientation_Second));
                nOrientation_Second(j) = length(ListOrientation_Second); 
            end
    end

    nOrientation_Second2 = nOrientation_Second./sampleRate;
    BinsNbSpikes_Second = nSpikesOrientation_Second./nOrientation_Second2;

    %  third quarter
    for j = 1:dBins;
        k = 0:dBins;
        ListOrientation_Third = find(Third_degf >= k(j)*6 & Third_degf < (k(j)+1)*6);
            if length(ListOrientation_Third) < 1;
                nOrientation_Third(j) = 1; 
                nSpikesOrientation_Third(j) = 0;
            else
                nSpikesOrientation_Third(j) = sum(Third_spksf(ListOrientation_Third));
                nOrientation_Third(j) = length(ListOrientation_Third); 
            end
    end

    nOrientation_Third2 = nOrientation_Third./sampleRate;
    BinsNbSpikes_Third = nSpikesOrientation_Third./nOrientation_Third2;

    %  fourth quarter
    for j = 1:dBins;
        k = 0:dBins;
        ListOrientation_Fourth = find(Fourth_degf >= k(j)*6 & Fourth_degf < (k(j)+1)*6);
            if length(ListOrientation_Fourth) < 1;
                nOrientation_Fourth(j) = 1; 
                nSpikesOrientation_Fourth(j) = 0;
            else
                nSpikesOrientation_Fourth(j) = sum(Fourth_spksf(ListOrientation_Fourth));
                nOrientation_Fourth(j) = length(ListOrientation_Fourth); 
            end
    end

    nOrientation_Fourth2 = nOrientation_Fourth./sampleRate;
    BinsNbSpikes_Fourth = nSpikesOrientation_Fourth./nOrientation_Fourth2;

    FirstSecond_corr = corrcoef(BinsNbSpikes_First',BinsNbSpikes_Second'); 
    FirstThird_corr = corrcoef(BinsNbSpikes_First',BinsNbSpikes_Third'); 
    FirstFourth_corr = corrcoef(BinsNbSpikes_First',BinsNbSpikes_Fourth'); 
    SecondThird_corr = corrcoef(BinsNbSpikes_Second',BinsNbSpikes_Third'); 
    SecondFourth_corr = corrcoef(BinsNbSpikes_Second',BinsNbSpikes_Fourth');
    ThirdFourth_corr = corrcoef(BinsNbSpikes_Third',BinsNbSpikes_Fourth');

    FourQuart_MeanCorr = (FirstSecond_corr(1,2) + FirstThird_corr(1,2) + FirstFourth_corr(1,2) + SecondThird_corr(1,2) + SecondFourth_corr(1,2) + ThirdFourth_corr(1,2))/6;
    
    % firing rate x HD polar plot for the data above
    polarplot = polar(BinsAngle([1:60 1]),BinsNbSpikes([1:60 1]),'b');
    set(polarplot, 'linewidth',3,'color','k');
    axis tight
%     set(gca,'XLim',[-16 16],'YLim',[-16 16]); %Set limits on polar axis to normalize across within-subject design
    set(0,'Showhiddenhandles','on')
    extrastuff = setdiff(get(gca,'children'),polarplot);
    delete(extrastuff)
    horizontal=line([-100 100],[0 0]);
    vertical=line([0 0],[-100 100]);
    set(horizontal,'linewidth',2,'color','k');
    set(vertical,'linewidth',2,'color','k');
    [filepath,filename] = fileparts(ReadData{i});
    saveas(polarplot,[filepath filesep filename '_polar.jpg']);
        
    keep('path', 'ReadData', 'sampleRate', 'xmin', 'xmax', 'ymin', 'ymax', 'minLED', 'maxLED', 'dBins', 'i','mean_vector_length', 'peak_Firing_Rate', 'preferred_Direction', 'overall_rate', 'Direct_infoContent', 'FirstSec_Corr', 'FourQuart_MeanCorr', 'Directional_Range_HalfWidth', 'Directional_Range_HalfWidth_bins');
    [filepath,filename] = fileparts(ReadData{i});
    save([filepath filesep filename '_HD_Properties.mat']);    
    
    close all    
    keep('path', 'ReadData');
end

clear all