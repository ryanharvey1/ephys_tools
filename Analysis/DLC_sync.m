%% DLC Tracking sync with neuralynx ephys

%% Synchronize behavior with electrophysiology
% Import DLC output
file = dir('*_dlc_tracking.csv');
dlc_res = readtable(file.name,'HeaderLines',2);

% Import Neuralynx video timestamps
ts = Nlx2MatVT([pwd,filesep,'VT1.nvt'],[1,0,0,0,0,0],0,1);

% Set timestamps to start of video 
if size(ts',1) ~= size(dlc_res,1)
    % size of each matrix
    size_mat =  [size(ts',1),size(dlc_res,1)];
    [~,idx]  =  min(size_mat);
    if idx == 1 % NVT is shorter 
        
    else % DLC is shorter (usually the case)
        timestamps = dlc_res.coords;
        red_coords = [dlc_res.x,dlc_res.y,dlc_res.likelihood];
        green_coords = [dlc_res.x_1,dlc_res.y_1,dlc_res.likelihood_1];
    
    
  nan(abs(diff(size_mat)),1)
    end
    
else
    
end


figure;
plot(frames.coords,frames.x,'k')
hold on;
idx = frames.likelihood < 0.99;
plot(frames.coords(idx),frames.x(idx),'*r')
