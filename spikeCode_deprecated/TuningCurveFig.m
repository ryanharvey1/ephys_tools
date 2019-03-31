
%Purpose of this code is to create a plot of all tuning curves for directionally modulated
%cells between stability sessions. 
%by L.Berkowitz March 2017

addpath(genpath('D:\ClarkP30_Recordings'));
addpath(genpath('F:\Users\reharvey\Place_Cell_Data\PAE_Rat\'));


CompileFilteredMatFiles

%plot 
figure (i+2), plot(BinsAngle3,smooth_BinsNbSpikes,'LineWidth',2,'color','k') %plot(BinsAngle3,BinsNbSpikes,'LineWidth',2,'color','k')
            axis tight;
            hold on
            xlim ([0 360]);
            box off
            fig3=figure (i+2);