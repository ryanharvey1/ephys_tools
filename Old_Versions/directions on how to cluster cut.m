%% Directions on how to cluster cut
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
hcList = {1 2 3 4 5 6 7 8};
preprocessTTUofL('F:\Users\BClarkLab\Desktop\2016-05-14_10-21-40','F:\Users\BClarkLab\Desktop\2016-05-14_10-21-40p',[], hcList);  

%% Check for Bad Channels
check_channels_kk('F:\Users\BClarkLab\Desktop\2016-05-14_10-21-40')
    % You will be presented with a series of graphs showing the activity on each channel of each tetrode.
    % Verify that all bad channels are checked, check or uncheck as necessary.
    % You will also need to indicate the number of cells on each tetrode using the dropdown menus over each tetrode. 
    % Estimated number of cells should have been determined during the experimental recording session.
    % When done selecting the number of units, click "Write File" and then "exit".
    
%% Run the Batch processing program.    
cd('F:\Users\BClarkLab\Desktop\2016-05-14_10-21-40p\TT')
RunClustBatch('batch.txt');

 %%
MClust % not working / Still no clu file
    %When cutting, cut only files found in the TT directory.

