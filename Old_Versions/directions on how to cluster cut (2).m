%% Directions on how to cluster cut

% From Preparing data for cluster cutting_may_06_2010 Ver3.pdf 
% set your current directory 
    %example directory
    % C:\Users\Ben Clark Lab\Documents\MATLAB\2016-05-12_15-15-14   
    % could also just use the file directory above (point and click)
edit preprocessTTUofL  
    % check to see if your correct template_filter is selcted
    % line 269 - 274
    % version 3.2 might have updated this ***********
    
    %also check to see if your largest .ntt is larger than 90mb (line 279)
%% running preprocess     
hcList = {1, 2, 3, 4, 5, 6, 7, 8};
preprocessTTUofL('F:\Users\BClarkLab\Desktop\2016-05-14_10-21-40','F:\Users\BClarkLab\Desktop\2016-05-14_10-21-40pp',[],hcList);  
    % are we really using hipp template?
    % Missing clu, klg, and model
%% Check for Bad Channels
check_channels_kk('F:\Users\BClarkLab\Desktop\2016-05-14_10-21-40')
    % You will be presented with a series of graphs showing the activity on each channel of each tetrode.
    % Verify that all bad channels are checked, check or uncheck as necessary.
    % You will also need to indicate the number of cells on each tetrode using the dropdown menus over each tetrode. 
    % Estimated number of cells should have been determined during the experimental recording session.
    % When done selecting the number of units, click "Write File" and then "exit".
    
%% Run the Batch processing program.    
cd('F:\Users\BClarkLab\Desktop\2016-05-14_10-21-40pp\TT')
RunClustBatch('batch.txt');

 %%
MClust % not working / Still no clu file
% Click "run klustakwik"
    % click "create/load FD..."
    % select NTT file
    % export clusters in the new window
    % unclick "run klustakwik"
    % click "create/load FD..." again
    % Select NTT file
    % another window will pop up. Click on cluster file you just created
    
    %  When cutting, cut only files found in the TT directory.
    % odd pop up window after you load your ntt file
    
    
%% Pipeline Ben and Aaron came up with
MClustPath = 'F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\MClust3.1UofL_c';
addpath(genpath(MClustPath));

input = 'F:\Users\BClarkLab\Desktop\2016-06-27_15-47-44';

output = 'F:\Users\BClarkLab\Desktop\2016-06-27_15-47-44p';

preprocessTTUofL(input,output);
%% 
check_channels_kk(output)
%%
cd 'F:\Users\BClarkLab\Desktop\2016-06-27_15-47-44p\TT'

RunClustBatch('Batch_KKwik.txt');

%%
MClust
% Click "run klustakwik"
    % click "create/load FD..."
    % select NTT file
    % export clusters in the new window
    % unclick "run klustakwik"
    % click "create/load FD..." again
    % Select NTT file
    % another window will pop up. Click on cluster file you just created
    
    %  When cutting, cut only files found in the TT directory.
    % odd pop up window after you load your ntt file
    %%
    beep on
    beep 
    beep