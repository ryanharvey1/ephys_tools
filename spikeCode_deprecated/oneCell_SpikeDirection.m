% SpikesDirection creates summary of data in form of spikes on path and
% tuning curve
%
% Input: 
%       - raw data output from Labview
%
% Output:
%       - figure showing spikes on path
%       - figure showing plot of HD cell tuning function
%       - figure showing plot of spikes x HD x timestamps
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
% Created by Ben C 2013; modified by Ben C Jan 2014

% identify path to data files
path = '/Users/bjclark/Desktop/Dropbox/EC_analysis/__RawData_&_Summary/__control_results_files_PaS/';
ReadData = FindFiles('bc217_s03_TT1_u2.txt', 'StartingDirectory', path);

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

% load xy data
data = importdata(ReadData{1});

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

% number of directional bins
dBins = 60;

% Binning of directional data can also be done using histc function
for i = 1:dBins;
    % number of instances where the head of the rat was found in each of the 60 possible orientations
    ListOrientation = find(anglef >= ((((2*pi)/dBins)/2)+(i-1)*2*pi/dBins) & anglef < ((((2*pi/dBins)/2)+(i)*2*pi/dBins)));
        if length(ListOrientation) < 1; % if the number is less than 1 (i.e., 0), it is hard set to 1 to avoid division by 0
           nOrientation(i) = 1; 
           nSpikesOrientation(i) = 0; % 0 is assined to the number of spikes for this orientation
        else
           nSpikesOrientation(i) = sum(spksf(ListOrientation));
           nOrientation(i) = length(ListOrientation);
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
ij= find(isinf(log_IC)); % find -Inf's (log(0)) and replace with 0's (based on code from McN lab)
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

for i = 1:dBins;
    ListOrientation_First = find(First_anglef >= ((((2*pi)/dBins)/2)+(i-1)*2*pi/dBins) & First_anglef < ((((2*pi/dBins)/2)+(i)*2*pi/dBins)));
        if length(ListOrientation_First) < 1;
           nOrientation_First(i) = 1; 
           nSpikesOrientation_First(i) = 0;
        else
           nSpikesOrientation_First(i) = sum(First_spksf(ListOrientation_First));
           nOrientation_First(i) = length(ListOrientation_First); 
    end
end

nOrientation_First2 = nOrientation_First./sampleRate;
BinsNbSpikes_First = nSpikesOrientation_First./nOrientation_First2;

for i = 1:dBins;
    ListOrientation_Sec = find(Sec_anglef >= ((((2*pi)/dBins)/2)+(i-1)*2*pi/dBins) & Sec_anglef < ((((2*pi/dBins)/2)+(i)*2*pi/dBins)));
        if length(ListOrientation_Sec) < 1;
           nOrientation_Sec(i) = 1; 
           nSpikesOrientation_Sec(i) = 0;
        else
           nSpikesOrientation_Sec(i) = sum(Sec_spksf(ListOrientation_Sec));
           nOrientation_Sec(i) = length(ListOrientation_Sec); 
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
for i = 1:dBins;
    ListOrientation_F = find(F_anglef >= ((((2*pi)/dBins)/2)+(i-1)*2*pi/dBins) & F_anglef < ((((2*pi/dBins)/2)+(i)*2*pi/dBins)));
        if length(ListOrientation_F) < 1;
           nOrientation_F(i) = 1; 
           nSpikesOrientation_F(i) = 0;
        else
           nSpikesOrientation_F(i) = sum(F_spksf(ListOrientation_F));
           nOrientation_F(i) = length(ListOrientation_F); 
    end
end

nOrientation_F2 = nOrientation_F./sampleRate;
BinsNbSpikes_F = nSpikesOrientation_F./nOrientation_F2;

%  second quarter
for i = 1:dBins;
    ListOrientation_S = find(S_anglef >= ((((2*pi)/dBins)/2)+(i-1)*2*pi/dBins) & S_anglef < ((((2*pi/dBins)/2)+(i)*2*pi/dBins)));
        if length(ListOrientation_S) < 1;
           nOrientation_S(i) = 1; 
           nSpikesOrientation_S(i) = 0;
        else
           nSpikesOrientation_S(i) = sum(S_spksf(ListOrientation_S));
           nOrientation_S(i) = length(ListOrientation_S); 
    end
end

nOrientation_S2 = nOrientation_S./sampleRate;
BinsNbSpikes_S = nSpikesOrientation_S./nOrientation_S2;

%  third quarter
for i = 1:dBins;
    ListOrientation_T = find(Third_anglef >= ((((2*pi)/dBins)/2)+(i-1)*2*pi/dBins) & Third_anglef < ((((2*pi/dBins)/2)+(i)*2*pi/dBins)));
        if length(ListOrientation_T) < 1;
           nOrientation_T(i) = 1; 
           nSpikesOrientation_T(i) = 0;
        else
           nSpikesOrientation_T(i) = sum(Third_spksf(ListOrientation_T));
           nOrientation_T(i) = length(ListOrientation_T); 
    end
end

nOrientation_T2 = nOrientation_T./sampleRate;
BinsNbSpikes_T = nSpikesOrientation_T./nOrientation_T2;

%  fourth quarter
for i = 1:dBins;
    ListOrientation_Fourth = find(Fourth_anglef >= ((((2*pi)/dBins)/2)+(i-1)*2*pi/dBins) & Fourth_anglef < ((((2*pi/dBins)/2)+(i)*2*pi/dBins)));
        if length(ListOrientation_Fourth) < 1;
           nOrientation_Fourth(i) = 1; 
           nSpikesOrientation_Fourth(i) = 0;
        else
           nSpikesOrientation_Fourth(i) = sum(Fourth_spksf(ListOrientation_Fourth));
           nOrientation_Fourth(i) = length(ListOrientation_Fourth); 
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

% % plot smoothed spike on path plot using red LED
figure (1), plot(rxs,rys, 'LineWidth',1,'color',[0,0,0]+0.8); 
set(gca,'YDir','reverse');
axis square tight
hold on
scatter(spks_xs, spks_ys, 45, 'filled', 'k');
box off
title('Spike (black dots) on Path (gray)');

% % plot firing rate x HD
% figure (2), subplot(3,1,1), plot(BinsAngle3,BinsNbSpikes,'k--o','LineWidth',1)
% axis tight;
% hold on
% xlim ([0 360]);
% box off
% title('Entire Session');
% ylabel('Firing Rate (spikes/sec)');
% 
% % plot firing rate x HD for first half of session
% figure (2), subplot(3,1,2), plot(BinsAngle3,BinsNbSpikes_First,'k--o','LineWidth',1)
% axis tight;
% hold on
% xlim ([0 360]);
% box off
% title('First Half of Session');
% ylabel('Firing Rate (spikes/sec)');
% 
% % plot firing rate x HD for second half of session
% figure (2), subplot(3,1,3), plot(BinsAngle3,BinsNbSpikes_Sec,'k--o','LineWidth',1)
% axis tight;
% hold on
% xlim ([0 360]);
% box off
% title('Second Half of Session');
% ylabel('Firing Rate (spikes/sec)');
% xlabel('Head Direction (degrees)');

% firing rate x HD polar plot for the data above
% polarplot = polar(BinsAngle,BinsNbSpikes,'b'), title('Entire Session');
figure (2), polarplot = polar(BinsAngle([1:end 1]),BinsNbSpikes([1:end 1]),'b'), title('Entire Session');
set(polarplot, 'linewidth',3,'color','k');
axis tight
set(0,'Showhiddenhandles','on')
extrastuff = setdiff(get(gca,'children'),polarplot);
delete(extrastuff)
horizontal=line([-100 100],[0 0]);
vertical=line([0 0],[-100 100]);
set(horizontal,'linewidth',2,'color','k');
set(vertical,'linewidth',2,'color','k');
PFR = uicontrol('Style','text','position',[585 544 90 45]);
set(PFR,'String',num2str(peak_Firing_Rate),'background','w','fontsize',34);
% [filepath,filename] = fileparts(ReadData{1});
% saveas(polarplot,[filepath filesep filename 'polar.eps']);
% figure (3), subplot(3,1,2), polar(BinsAngle,BinsNbSpikes_First,'b'), title('First Half of Session');
% figure (3), subplot(3,1,3), polar(BinsAngle,BinsNbSpikes_Sec,'b'), title('Second Half of Session');
% figure (4), subplot(4,1,1), polar(BinsAngle,BinsNbSpikes_F,'b'), title('First Quarter of Session');
% figure (4), subplot(4,1,2), polar(BinsAngle,BinsNbSpikes_S,'b'), title('Second Quarter of Session');
% figure (4), subplot(4,1,3), polar(BinsAngle,BinsNbSpikes_T,'b'), title('Third Quarter of Session');
% figure (4), subplot(4,1,4), polar(BinsAngle,BinsNbSpikes_Fourth,'b'), title('Fourth Quarter of Session');

% % spikes x direction x timestamp
% figure (5), plot(degHeading, 'LineWidth',1,'color',[0,0,0]+0.8); 
% axis tight
% hold on
% ylim ([0 360]);
% scatter(tSpks, dSpks, 35, 'filled', 'b');
% box off
% title('Spikes (blue dots); Head Direction (gray lines)');
% xlabel('Timestamps');

% clear all variables bu measuresof HD
keep('mean_vector_length', 'peak_Firing_Rate', 'preferred_Direction', 'overall_rate', 'Direct_infoContent', 'FirstSec_Corr', 'FourQuart_MeanCorr', 'Directional_Range_HalfWidth', 'Directional_Range_HalfWidth_bins', 'BinsNbSpikes');

