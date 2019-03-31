path = '/Users/bjclark/Desktop/Dropbox/place cell analysis/Raw Data/RH026 10-28-2013/Spike Timestamps';
ReadData = FindFiles('TT2- RH_SS_09.txt', 'StartingDirectory', path);

    %load spike train data from .ts file outputs (.ts.r)
    ts = dlmread(ReadData{1},'\t',0,0); %load spike data
    
    %pad ts to have enough cells to run iterative code
    pad=zeros(150,1);
    ts_pad=[ts;pad];
    
    %convert spikes from 10s of usec to msec
    ts_msec = ts_pad./1000;
    
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
        else ISIx >= 20;
        end
    end
    
    %remove all negative values, zeros, and values over 500
    ISIfinal0=ISI_list;
    ISIfinal = ISIfinal0(ISIfinal0~=0);
    ISIfinal(ISIfinal <= 0) = [];
    ISIfinal(ISIfinal > 20) = [];
    
    %bin ISI by delay
    [ISIbins]=hist(ISIfinal,(1:1000));
%         
%     %smooth autocorrelation curve and normalize to 1
%     ISIsmooth = runline(ISIbins, 10, 1);
%     ISIsmooth(ISIsmooth < 0) = [0];
%     maxISIsmooth = max(ISIsmooth);
%     normISIsmooth = ISIsmooth./maxISIsmooth;
    
%     %calculate a theta ratio (first trough/second peak)
%     minISI = min(normISIsmooth(25:105));
%     maxISI = max(normISIsmooth(95:175));
%     ThetaRatio = (maxISI-minISI);
    
%     %fit a sinusoid to the autocorrelation
%     xSinu=[0.00628318:0.00628318:3.14159]';
%     opts = fitoptions('Method', 'NonlinearLeastSquares');
%     opts.Lower = [-Inf -Inf -Inf 4];
%     opts.StartPoint = [0 0 0 7];
%     opts.Upper = [Inf Inf Inf 12];
%     [f, G]=fit(xSinu,normISIsmooth,'fourier1',opts);
%     Frequency=f.w;
%     Fit=G.rsquare;
    
    %plot the autocorrelation
    fig1 = area(normISIsmooth,'facecolor','k'); %smoothed and normalized autocorrelation
    axis square tight
    set(gca,'XTick',[0 5 10 15 20],'YTick',[0 0.5 1.0]);
    box off
    
%     %save data and fig
%     [filepath, filename] = fileparts(ReadData{i});
%     saveas(fig1,[filepath filesep filename 'ISItheta.jpg']);
%     save([filepath filesep filename '_ISITheta.mat']);
%     close all
