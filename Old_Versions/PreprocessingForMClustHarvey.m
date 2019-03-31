%% PreprocessingForMClustHarvey 
% This first section will run through a parent folder and preprocess all
% folders(sessions)inside so that they will be ready for MClust
% Ryan H June 2016
% Wishlist: 1.When directed at a folder, only process new sessions. 2.Convert to function to protect code changes.

clear
addpath(genpath('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\MClust3.1UofL_c'));

disp('Locate Parent Directory');
workingDir= uigetdir; % opens up file window
fldr_list = dir(workingDir); 

for i = 1:length(fldr_list);
    if fldr_list(i).isdir && fldr_list(i).name(1) ~= '.' %&& any(regexp(fldr_list(i).name,'p$')) && any(regexp(fldr_list(i+1).name,'p$'))==0 ||==[]% removes '.' & don't run p files or files already ran
        cd([workingDir filesep fldr_list(i).name]);
        input_ = [workingDir filesep fldr_list(i).name];
        output =[input_,'p'];
        try
            preprocessTTUofL(input_,output);
        catch
        end
    end
end
disp(' ');
beep on
beep
% disp('You have 15 seconds to answer next question');
% disp(' ');
% disp('Do you want to cycle through files automatically? y / n:');
% choice = input(' ', 's');
% pause(15); % need a pause in order to write and save file to correct pwd in GUI window

% if choice=='y'; % automatically go through each GUI window
%

% FileInfo = dir('D:\Place_Cell_Data\PAE_Rat\Test_83016\2016-08-06_12-34-42p\TT\TT01c.ntt');
% [Y, M, D, H, MN, S] = datevec(FileInfo.datenum);
% c=clock;
% if Y==c(1,1) && M==c(1,2) && D==c(1,3) && H==c(1,4) 
%     disp('cool!')
% end
% 
% %
    pfolders= FindFiles('*p','StartingDirectory', workingDir);
    for ii=1:length(pfolders);
%         FileInfo = dir([pfolders(ii),'\TT\TT01c.ntt']);
%         [Y, M, D, H, MN, S] = datevec(FileInfo.datenum);
%         c=clock;
%         if Y==c(1,1) && M==c(1,2) && D==c(1,3) && H==c(1,4) 
            try
                check_channels_kk(pfolders{ii});
            catch
            end
%         end
    pause(15);
    end
  
% if choice=='n'; % go through each GUI window manually (Not tested)
%     for iiii=1:length(FindFiles('*p','StartingDirectory', workingDir));
%     pfolders= uigetdir;
%         try
%             check_channels_kk(pfolders);
%         catch
%         end
%     end
% end
% end
pfolders= FindFiles('*p','StartingDirectory', workingDir);
    for iii=1:length(pfolders);
        try
            fullname=fullfile(pfolders{iii},'TT');
            cd(fullname);
            RunClustBatch('Batch_KKwik.txt');
        catch
        end
    end
% cd('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory');
% Emailcode();


%% Manual Session by session preprocessing ***USE THIS IN MOST CASES***
clear
addpath(genpath('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\MClust3.1UofL_c'));
input_= uigetdir; 
output_ =[input_,'p'];
preprocessTTUofL(input_,output_,'hipp',(1:8));

check_channels_kk(output_);
pause(10); % pause to write and okay saved batch file

cd([output_ filesep 'TT']);
RunClustBatch('Batch_KKwik.txt');

  
%% opens MClust
clear
addpath(genpath('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\MClust3.1UofL_c'));
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
    
% Cutting/merging with MClust: Instructions, search ('Preparing Data for Cluster Cutting')



%% Post-Processing MClust Output into 1 figure with all spikes over time
NSMA_Toolbox = ('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\NSMA_Toolbox');
    addpath(genpath(NSMA_Toolbox));
% addfolders = 'F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\Cell analysis\Nlx2MatVt';
%     addpath(genpath(addfolders)); 
cd('D:\Place_Cell_Data\PAE_Rat\Test_83016\2016-08-25_16-13-19p\TT');
    ls '*.t'
tfile=FindFiles('*.t');
% tfile = ('D:\Place_Cell_Data\PAE_Rat\RH16\2016-07-02_18-35-24p\TT\TT01c_1__');
S =LoadSpikes(tfile);
    whos s
        % 
        %     cd('D:\Place_Cell_Data\PAE_Rat\RH16\2016-07-02_18-35-24p')
        %         VTfile = FindFiles('*.nvt');
        %         rmpath(genpath(NSMA_Toolbox));
        %     [VTTimeStamps, ExtractedX, ExtractedY, ExtractedAngle] = Nlx2MatVT(VTfile{1}, [1 1 1 1 0 0], 0, 1); 
        %         NSMA_Toolbox = ('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\NSMA_Toolbox');
        %         addpath(genpath(NSMA_Toolbox));
        %     [startrun, stoprun, Vel_cmps, Vel_phi] = FindRatRunsXY(ExtractedX, ExtractedY, 2.1, 5);
        %     X_run = Restrict(X, startrun, stoprun);
        %     nCells = length(S);
        %         for c=1:nCells
        %             S_run{c} = Restrict(S{c}, startrun, stoprun);
        %         end
        %     x_spikes{c} = interp1(Range(X_run), Data(X_run), Range(S_run{c}), 'linear');
        %     y_spikes{c} = interp1(Range(Y_run), Data(Y_run), Range(S_run{c}), 'linear');
        %     plot(x_spikes{c}, y_spikes{c}, 'r.')
Q = MakeQfromS(S,1000);
    whos Q
QD = Data(Q);
    whos QD
QD = QD';
    whos QD
figure; imagesc(QD)





