% spikes_shuffle shifts the spike series relative to the time series
%
% Input: 
%       - raw data output from Labview
%
% Output:
%       - rSHUFF = mean vector length for shuffled spike series
%       - corrSHUFF
%
% Created by Ben C 2014

% identify path to data files
path = '/Users/bjclark/Documents/MATLAB/test_data/spike_data/real_tests/2010-08-05_18-17-44/';
ReadData = FindFiles('bc219*.txt', 'StartingDirectory', path);

% video frame rate in Hz
sampleRate = 60;

% maze parameters
xmin = 0;
xmax = 255;
ymin = 0;
ymax = 255;

% min and max LED distance parameters
minLED = 4;
maxLED = 39;

% number of directional bins
dBins = 60;
shuff_interv = dBins*20;

for i = 1:length(ReadData);
    
    data = importdata(ReadData{i}); % load xy data

    % extract coords, spikes, angle, and direction from ReadData output
    rx = data.data(:,2); % red x-coords
    ry = data.data(:,3); % red y-coords
    gx = data.data(:,4); % green x-coords
    gy = data.data(:,5); % green y-coords
    spks = data.data(:,6); % spikes
    angle = data.data(:,10); % angle
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

    % extract spike series and angle
    spksf = rgmmFILT(:,5);
    anglef = rgmmFILT(:,6);

    % xcorr of tuning function between four quarters of session
    sessLength_samples = numel(anglef);
    F_End = round(sessLength_samples/4);
    S_Start = F_End + 1;
    S_End = F_End*2;
    Third_Start = S_End + 1;
    Third_End = F_End*3;
    Fourth_Start = Third_End + 1;
    Fourth_End = sessLength_samples;

    F_anglef = anglef(1:F_End);
    S_anglef = anglef(S_Start:S_End);
    Third_anglef = anglef(Third_Start:Third_End);
    Fourth_anglef = anglef(Fourth_Start:Fourth_End);

    rSHUFF = [];
    corrSHUFF = [];
        for x = 1:400;
            spkSHUFFi = circshift(spksf,randi([-shuff_interv shuff_interv],1));
            for j = 1:dBins;
                ListOrientation = find(anglef >= ((((2*pi)/dBins)/2)+(j-1)*2*pi/dBins) & anglef < ((((2*pi/dBins)/2)+(j)*2*pi/dBins)));
                    if length(ListOrientation) < 1; 
                        nOrientation(j) = 1; 
                        nSpikesOrientation(j) = 0; % 0 is assined to the number of spikes for this orientation
                    else
                        nSpikesOrientation(j) = sum(spkSHUFFi(ListOrientation));
                        nOrientation(j) = length(ListOrientation); 
                    end
            end
            nOrientation2 = nOrientation./sampleRate;
            BinsNbSpikes = nSpikesOrientation./nOrientation2;
            BinsAngle = (0+((2*pi/(dBins)/2)):2*pi/(dBins):(2*pi)-((2*pi/(dBins)/2)));
            spksi = find(spkSHUFFi == 1);
            aSpks = anglef(spksi);
            r = circ_r(BinsAngle',BinsNbSpikes',circ_ang2rad(6));
            rSHUFF = [rSHUFF; r];
    
            F_spksf = spkSHUFFi(1:F_End);
            S_spksf = spkSHUFFi(S_Start:S_End);
            Third_spksf = spkSHUFFi(Third_Start:Third_End);
            Fourth_spksf = spkSHUFFi(Fourth_Start:Fourth_End);
    
            % first quarter
            for k = 1:dBins;
                ListOrientation_F = find(F_anglef >= ((((2*pi)/dBins)/2)+(k-1)*2*pi/dBins) & F_anglef < ((((2*pi/dBins)/2)+(k)*2*pi/dBins)));
                    if length(ListOrientation_F) < 1;
                        nOrientation_F(k) = 1; 
                        nSpikesOrientation_F(k) = 0;
                    else
                        nSpikesOrientation_F(k) = sum(F_spksf(ListOrientation_F));
                        nOrientation_F(k) = length(ListOrientation_F); 
                    end
            end

            nOrientation_F2 = nOrientation_F./sampleRate;
            BinsNbSpikes_F = nSpikesOrientation_F./nOrientation_F2;

            %  second quarter
            for l = 1:dBins;
                ListOrientation_S = find(S_anglef >= ((((2*pi)/dBins)/2)+(l-1)*2*pi/dBins) & S_anglef < ((((2*pi/dBins)/2)+(l)*2*pi/dBins)));
                    if length(ListOrientation_S) < 1;
                        nOrientation_S(l) = 1; 
                        nSpikesOrientation_S(l) = 0;
                    else
                        nSpikesOrientation_S(l) = sum(S_spksf(ListOrientation_S));
                        nOrientation_S(l) = length(ListOrientation_S); 
                    end
            end

            nOrientation_S2 = nOrientation_S./sampleRate;
            BinsNbSpikes_S = nSpikesOrientation_S./nOrientation_S2;

            %  third quarter
            for m = 1:dBins;
                ListOrientation_T = find(Third_anglef >= ((((2*pi)/dBins)/2)+(m-1)*2*pi/dBins) & Third_anglef < ((((2*pi/dBins)/2)+(m)*2*pi/dBins)));
                    if length(ListOrientation_T) < 1;
                        nOrientation_T(m) = 1; 
                        nSpikesOrientation_T(m) = 0;
                    else
                        nSpikesOrientation_T(m) = sum(Third_spksf(ListOrientation_T));
                        nOrientation_T(m) = length(ListOrientation_T); 
                    end
            end

            nOrientation_T2 = nOrientation_T./sampleRate;
            BinsNbSpikes_T = nSpikesOrientation_T./nOrientation_T2;

            %  fourth quarter
            for n = 1:dBins;
                ListOrientation_Fourth = find(Fourth_anglef >= ((((2*pi)/dBins)/2)+(n-1)*2*pi/dBins) & Fourth_anglef < ((((2*pi/dBins)/2)+(n)*2*pi/dBins)));
                    if length(ListOrientation_Fourth) < 1;
                        nOrientation_Fourth(n) = 1; 
                        nSpikesOrientation_Fourth(n) = 0;
                    else
                        nSpikesOrientation_Fourth(n) = sum(Fourth_spksf(ListOrientation_Fourth));
                        nOrientation_Fourth(n) = length(ListOrientation_Fourth); 
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
            corrSHUFF = [corrSHUFF; FourQuart_MeanCorr]; 
        end
        
%         clear all non essential variables
        keep('rSHUFF', 'corrSHUFF', 'path', 'ReadData', 'sampleRate', 'xmin', 'xmax', 'ymin', 'ymax', 'minLED', 'maxLED', 'dBins', 'shuff_interv', 'i');
        [filepath, filename] = fileparts(ReadData{i});
        save([filepath filesep filename '_SHUFF.mat']);  
end


% % shuffled ditribution for firing rate
% percNinetyFive_MVR = prctile(rSHUFF,95);
% percNinetyNine_MVR = prctile(rSHUFF,99);
% subplot(2,1,1), hist(rSHUFF,100);
% hold on
% line([percNinetyNine_MVR percNinetyNine_MVR], [0 30], 'LineWidth',1, 'color', 'r')
% line([percNinetyFive_MVR percNinetyFive_MVR], [0 30], 'LineWidth',1, 'color', 'b')
% title('Shuffled Distribution of R-Values Calculated from Binned Firing Rate (red = 99th, blue = 95th)');
% 
% % shuffled ditribution for corr
% percNinetyFive_CORR = prctile(corrSHUFF,95);
% percNinetyNine_CORR = prctile(corrSHUFF,99);
% subplot(2,1,2), hist(corrSHUFF,100);
% hold on
% line([percNinetyNine_CORR percNinetyNine_CORR], [0 30], 'LineWidth',1, 'color', 'r')
% line([percNinetyFive_CORR percNinetyFive_CORR], [0 30], 'LineWidth',1, 'color', 'b')
% title('Shuffled Distribution of Mean Corr (red = 99th, blue = 95th)');

% 
% path = '/Users/bjclark/Desktop/Dropbox/ec recordings/__RawData_&_Summary/ben_data/control_sham/BC219 - TT 12HD 4BC 1GC 2GCorPC/2010-08-05_18-17-44/';
% ReadData2 = FindFiles('bc219_s9_TT3_u1_c1.txt', 'StartingDirectory', path);


% clear ReadData2 path




