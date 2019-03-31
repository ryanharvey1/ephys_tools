% HD_theta_skip

% opens .ts_R text files containing HD cell data and plots autoCorrISI
%
%
% Ryan E Harvey 2018
%
% cd to data and get file names
% clear;clc;close all
com=which('HD_theta_skip');
com=strsplit(com,filesep);
basedir=[com{1},filesep,'Users',filesep,com{3},filesep,'GoogleDrive',filesep,'MatlabDir'];
addpath([basedir,filesep,'BClarkToolbox',filesep, 'Analysis']);

% if ismac
%     path='/Users/ryanharvey/Downloads/HeadDIrectionCells_LauraRyan';
% else
    path='D:\Projects\Multi_Region_HD\HeadDIrectionCells_wRfiles';
% end

% path='/Users/lauraberkowitz/Google Drive/Manuscripts/In Progress/Ben_HDCProject/Data';

cd(path)
% files=dir( '**/*.ts_R');
files=dir('**/*_Timestamps'); %files=dir('**/*_Timestamps*.*')
for i=1:length(files)
    folders{i,:}=[files(i).folder,filesep,files(i).name];
end

disp(folders)

for i=1:length(folders)
    area=strsplit(folders{i},filesep);
    area=area{end};
    area=erase(area,'_Timestamps');
    areas{i}=erase(area,'_');
end
disp(areas)
clear area files 


for a=1:length(folders)
    cd(folders{a})
    filenames=dir;
    filenames={filenames.name}';
    filenames(strcmp(filenames,'..') | strcmp(filenames,'.') )=[];
    disp([num2str(length(filenames)),' ',areas{a},' cells'])
    for i=1:length(filenames)
        % OPEN TEXT FILE AND CLEAN DATA
        disp(['Running:  ',areas{a},filenames{i}])
        
        fileID=fopen(filenames{i},'r');
        dataArray = textscan(fileID, '%f%*s%[^\n\r]', 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,2-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        fclose(fileID);
        ts=[dataArray{1:end-1}];
        clearvars fileID dataArray ans;
        
        ts(ts<0)=[];
        
        % spike times to seconds
        spk=(ts./100)./1000;
        
        % compute mle rhythmicity
        [params,confidence,stats,everything]=mle_rhythmicity(spk,max(spk),'max_lag',0.5,'t_axis',linspace(0,.5,51),'plotit',0);
        close
        mle(i,1)=params.a;
        %% theta skip

        max_lag=.5;
        t_bin=0.010; % 10ms
        if t_bin / mod(max_lag, t_bin) ~= 2 % set lags so it is 'even' (odd number of coefficients and zero centered')
            max_lag = t_bin*floor(max_lag/t_bin)+.5*t_bin;
        end
        [skipcor, lag] = CrossCorr(spk, 'lag', [-max_lag max_lag], 'binsize', t_bin, 'norm', 'prob');

        thetaskip(i,1)=ThetaSkipping(skipcor,lag,t_bin); % theta skip
        
        thetaindex(i,1)=thetamod(skipcor,lag); %theta index
        
        smoothCor=smoothdata(skipcor,'gaussian',5);
        
        coherence(i,1)=corr(skipcor,smoothCor); %coherence smoothed versus raw
        
        ids{i,1}=filenames{i};
        
        thetaInfo=[thetaskip(i,1),thetaindex(i,1), size(spk,1),coherence(i,1)];
        
%         fig=figure; fig.Color=[1 1 1];
%         plot(lag,skipcor,'LineWidth',2, 'color','k'); hold on;
%         axis tight
%         grid on
%         box off
%         plot(lag,smoothCor,'LineWidth',2,'color','r');
%         title(sprintf('Theta Skipping: %4.2f, Theta Index: %4.2f, Num Spikes: %d, Coherence: %4.2f',thetaInfo))
%         
% %         print(fig,'-dpng', '-r80',['C:\Users\ryanh\Dropbox\school work\UNM\Lab\Projects\Multi_Region_HDcells\Figures\Theta_Skipping\autocorr',filesep,areas{a},'_',ids{i,1},'.png'])
%         print(fig,'-dpng', '-r80',['/Users/lauraberkowitz/Google Drive/Manuscripts/In Progress/Ben_HDCProject/Data/thetaFigs',filesep,areas{a},'_',ids{i,1},'.png'])
% 
%         close
    end
    HDdata.(areas{a}).id=ids;
    
    HDdata.(areas{a}).thetaskip=thetaskip;
    HDdata.(areas{a}).thetaCoherence=coherence;
    HDdata.(areas{a}).mle=mle;
    HDdata.(areas{a}).thetaindex=thetaindex;


    clear  ids  thetaskip coherence thetaindex mle
    
end

%% PLOT

% theta skipping by layer
clear data
data.ATN=HDdata.ATN.thetaskip;
data.PoSDeep=HDdata.PoSDeep.thetaskip;
data.PoSSup=HDdata.PoSSup.thetaskip;
data.MECDeepLayers=HDdata.MECDeepLayers.thetaskip;
data.MECSupLayers=HDdata.MECSupLayers.thetaskip;
data.PaSDeep=HDdata.PaSDeep.thetaskip;
data.PaSSup=HDdata.PaSSup.thetaskip;

fig=figure;fig.Color=[1 1 1];
ECDF_plot(data,'Theta Skipping');


% mle by layer
clear data
data.ATN=HDdata.ATN.mle;
data.PoSDeep=HDdata.PoSDeep.mle;
data.PoSSup=HDdata.PoSSup.mle;
data.MECDeepLayers=HDdata.MECDeepLayers.mle;
data.MECSupLayers=HDdata.MECSupLayers.mle;
data.PaSDeep=HDdata.PaSDeep.mle;
data.PaSSup=HDdata.PaSSup.mle;

fig=figure;fig.Color=[1 1 1];
ECDF_plot(data,'mle');

% thetaindex by layer
clear data
data.ATN=HDdata.ATN.thetaindex;
data.PoSDeep=HDdata.PoSDeep.thetaindex;
data.PoSSup=HDdata.PoSSup.thetaindex;
data.MECDeepLayers=HDdata.MECDeepLayers.thetaindex;
data.MECSupLayers=HDdata.MECSupLayers.thetaindex;
data.PaSDeep=HDdata.PaSDeep.thetaindex;
data.PaSSup=HDdata.PaSSup.thetaindex;

fig=figure;fig.Color=[1 1 1];
ECDF_plot(data,'thetaindex');

% coherence by layer
clear data
data.ATN=HDdata.ATN.thetaCoherence;
data.PoSDeep=HDdata.PoSDeep.thetaCoherence;
data.PoSSup=HDdata.PoSSup.thetaCoherence;
data.MECDeepLayers=HDdata.MECDeepLayers.thetaCoherence;
data.MECSupLayers=HDdata.MECSupLayers.thetaCoherence;
data.PaSDeep=HDdata.PaSDeep.thetaCoherence;
data.PaSSup=HDdata.PaSSup.thetaCoherence;

fig=figure;fig.Color=[1 1 1];
ECDF_plot(data,'thetaCoherence');

%____________________________Local Function________________________________

function thetaindex=thetamod(cor,lag)

if length(lag)~=length(cor)
    lag=linspace(-.5,.5,length(cor));
end

cor = cor/max(cor(lag~=0))*2-1;
if any(cor>1)
   cor=rescale(cor,-1,1);
end
% The power spectrum of the temporal autocorrelograms was 
% assessed by computing the fast Fourier transform (FFT) of the 
% autocorrelogram, and calculating the square of the FFT magnitude
S = abs(fft(cor)).^2;
df = 1/range(lag); fNQ = 1/mode(diff(lag))/2;
f = 0:df:fNQ;
S = S(1:length(f));
% The power spectrum was smoothed with a 2-Hz rectangular window
S = smooth(S,2/df);

% and the peak value in the 4-12 Hz band was identified
peak = f(f>=4&f<=12);
[~,i] = max(S(f>=4&f<=12));
peak = peak(i);

% A neuron was defined as theta)modulated if the mean power within 1)Hz of 
% each side of the peak in the 5�11 Hz frequency range was at least 5 times
% greater than the mean spectral power between 0 Hz and 50 Hz
thetaindex = mean(S(abs(f-peak)<1))/mean(S(f<50));

end


