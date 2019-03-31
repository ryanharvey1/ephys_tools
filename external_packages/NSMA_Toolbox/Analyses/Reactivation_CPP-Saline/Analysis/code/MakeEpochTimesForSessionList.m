%% prints 3 formated strings to copy and paste into SessionList.m template
%% needs the epoch start-end times in workspace
%% e.g :  s.epochs.sleep1 = [14969140.37,20969140.37];


disp(['s.epochs.sleep1 = [' num2str(ts_start_s1) ',' num2str(ts_end_s1) '];'])
disp(['s.epochs.maze1  = [' num2str(ts_start_b1) ',' num2str(ts_end_b1) '];'])
disp(['s.epochs.sleep2 = [' num2str(ts_start_s2) ',' num2str(ts_end_s2) '];'])
