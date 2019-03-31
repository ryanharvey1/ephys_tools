function normISIsmooth=autoCorrISI(ts,fig)
% autoCorrISI: plots auto correlation 
%
% Input: ts
%           timestamps in 10s of secconds
%
% 
pad=zeros(150,1);
ts_pad=[ts;pad];

%convert spikes from 10s of usec to msec
ts_msec = ts_pad./100;

ISI_list=[];
ISIx=[];
ISIy=[];
y=([2:101])';
%subtract each ts from all subsequent ts to generate ISI list
for j = 1:(length(ts_pad)-149);
    ISIx = ts_msec(j+1)-ts_msec(j);
    ISI_list = [ISI_list; ISIx];
    if ISIx >= 0 && ISIx <= 501;
        ISIy = ts_msec(j+y)-ts_msec(j);
        ISI_list = [ISI_list; ISIy];
    else
        ISIx >= 500;
    end
end

%remove all negative values, zeros, and values over 500
ISIfinal0=ISI_list;
ISIfinal = ISIfinal0(ISIfinal0~=0);
ISIfinal(ISIfinal <= 0) = [];
ISIfinal(ISIfinal > 500) = [];

%bin ISI by delay
[ISIbins]=hist(ISIfinal,(1:500));

%smooth autocorrelation curve and normalize to 1
ISIsmooth = runline(ISIbins,10,1);
ISIsmooth(ISIsmooth < 0) = 0;
maxISIsmooth = max(ISIsmooth);
normISIsmooth = ISIsmooth./maxISIsmooth;

%plot the autocorrelation
if fig==1
    area(normISIsmooth,'facecolor','k'); %smoothed and normalized autocorrelation
    axis square tight
    set(gca,'XTick',[0 100 200 300 400 500],'YTick',[0 0.5 1.0]);
    box off
end
end