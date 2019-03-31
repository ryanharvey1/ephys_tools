function [ ThetaStats ] = ThetaPower( EEG_DownSampledData )
%ThetaPower Summary of this function goes here
% Input:
%         - Downsampled LFP
%         
% Output:
%         - theta power from 4-8Hz (raw from the spectrogram)
%         - theta ratio (theta power [4-12Hz]/delta power [1-4Hz])
%         - theta frequency (frequency between 4-8Hz with highest power)
%
% Created by Shawn Winter September 2014
% Adapted by Ryan E Harvey April 2017
%
[~,F,~,P] = spectrogram(EEG_DownSampledData,hanning(1000),512,0:0.25:20,1000);

% create variables for theta (6-10Hz) and delta (2-4Hz) power over time
theta = F(:,1) >= 4 & F(:,1) <= 12;
thetaP_all = P(theta,:);
delta = F(:,1) >= 1 & F(:,1) <= 4;
deltaP_all = P(delta,:);

%calculate theta ratio (theta power/delta power) over time 
thetaP_allm1 = mean(thetaP_all); %instantaneous power calculated from the mean
thetaP_allm2 = mean(thetaP_allm1); %mean overall power
deltaP_allm1 = mean(deltaP_all); %instantaneous power calculated from the mean
deltaP_allm2 = mean(deltaP_allm1); %mean overall power
ThetaRatio = thetaP_allm2/deltaP_allm2; %theta ratio from means

% Find the max power in each bin, eliminate all bins with max power outside of theta range
[MaxOverall_P,MaxOverall_L] = max(P);
MaxTheta_row = find(MaxOverall_L >= 17 & MaxOverall_L <= 49);
MaxTheta_P = MaxOverall_P(:,MaxTheta_row);
MaxTheta_L = MaxOverall_L(:,MaxTheta_row);

% Convert bin location to frequency and calculate overall mean theta frequency
MaxTheta_F=MaxTheta_L;
MaxTheta_F(MaxTheta_F==1)=2;
MaxTheta_F(MaxTheta_F==2)=2.25;
MaxTheta_F(MaxTheta_F==3)=2.5;
MaxTheta_F(MaxTheta_F==4)=2.75;
MaxTheta_F(MaxTheta_F==5)=3;
MaxTheta_F(MaxTheta_F==6)=3.25;
MaxTheta_F(MaxTheta_F==7)=3.5;
MaxTheta_F(MaxTheta_F==8)=3.75;
MaxTheta_F(MaxTheta_F==9)=4;
MaxTheta_F(MaxTheta_F==10)=4.25;
MaxTheta_F(MaxTheta_F==11)=4.5;
MaxTheta_F(MaxTheta_F==12)=4.75;
MaxTheta_F(MaxTheta_F==13)=5;
MaxTheta_F(MaxTheta_F==14)=5.25;
MaxTheta_F(MaxTheta_F==15)=5.5;
MaxTheta_F(MaxTheta_F==16)=5.75;
MaxTheta_F(MaxTheta_F==17)=6;
MaxTheta_F(MaxTheta_F==18)=6.25;
MaxTheta_F(MaxTheta_F==19)=6.5;
MaxTheta_F(MaxTheta_F==20)=6.75;
MaxTheta_F(MaxTheta_F==21)=7;
MaxTheta_F(MaxTheta_F==22)=7.25;
MaxTheta_F(MaxTheta_F==23)=7.5;
MaxTheta_F(MaxTheta_F==24)=7.75;
MaxTheta_F(MaxTheta_F==25)=8;
MaxTheta_F(MaxTheta_F==26)=8.25;
MaxTheta_F(MaxTheta_F==27)=8.5;
MaxTheta_F(MaxTheta_F==28)=8.75;
MaxTheta_F(MaxTheta_F==29)=9;
MaxTheta_F(MaxTheta_F==30)=9.25;
MaxTheta_F(MaxTheta_F==31)=9.5;
MaxTheta_F(MaxTheta_F==32)=9.75;
MaxTheta_F(MaxTheta_F==33)=10;
MaxTheta_F(MaxTheta_F==34)=10.25;
MaxTheta_F(MaxTheta_F==35)=10.5;
MaxTheta_F(MaxTheta_F==36)=10.75;
MaxTheta_F(MaxTheta_F==37)=11;
MaxTheta_F(MaxTheta_F==38)=11.25;
MaxTheta_F(MaxTheta_F==39)=11.5;
MaxTheta_F(MaxTheta_F==40)=11.75;
MaxTheta_F(MaxTheta_F==41)=12;
MaxTheta_F(MaxTheta_F==42)=12.25;
MaxTheta_F(MaxTheta_F==43)=12.5;
MaxTheta_F(MaxTheta_F==44)=12.75;
MaxTheta_F(MaxTheta_F==45)=13;
MaxTheta_F(MaxTheta_F==46)=13.25;
MaxTheta_F(MaxTheta_F==47)=13.5;
MaxTheta_F(MaxTheta_F==48)=13.75;
MaxTheta_F(MaxTheta_F==49)=14;
MaxTheta_F(MaxTheta_F==50)=14.25;
MaxTheta_F(MaxTheta_F==51)=14.5;
MaxTheta_F(MaxTheta_F==52)=14.75;
MaxTheta_F(MaxTheta_F==53)=15;
MaxTheta_F(MaxTheta_F==54)=15.25;
MaxTheta_F(MaxTheta_F==55)=15.5;
MaxTheta_F(MaxTheta_F==56)=15.75;
MaxTheta_F(MaxTheta_F==57)=16;
MaxTheta_F(MaxTheta_F==58)=16.25;
MaxTheta_F(MaxTheta_F==59)=16.5;
MaxTheta_F(MaxTheta_F==60)=16.75;
MaxTheta_F(MaxTheta_F==61)=17;
MaxTheta_F(MaxTheta_F==62)=17.25;
MaxTheta_F(MaxTheta_F==63)=17.5;
MaxTheta_F(MaxTheta_F==64)=17.75;
MaxTheta_F(MaxTheta_F==65)=18;
MaxTheta_F(MaxTheta_F==66)=18.25;
MaxTheta_F(MaxTheta_F==67)=18.5;
MaxTheta_F(MaxTheta_F==68)=18.75;
MaxTheta_F(MaxTheta_F==69)=19;
MaxTheta_F(MaxTheta_F==70)=19.25;
MaxTheta_F(MaxTheta_F==71)=19.5;
MaxTheta_F(MaxTheta_F==72)=19.75;
MaxTheta_F(MaxTheta_F==73)=20;
MeanThetaFreq = nanmean(MaxTheta_F);

% Put into Struct
% ThetaStats.MaxTheta_P=MaxTheta_P;
ThetaStats.MeanThetaFreq=MeanThetaFreq;
ThetaStats.MeanOverallPow=thetaP_allm2;
ThetaStats.ThetaRatio=ThetaRatio;
end

