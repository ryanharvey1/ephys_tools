function backup(BackUpSourceDirs,DestDir)
%file backup program
%Written by Matt Becker
%email:  mrmattbecker@hotmail.com
%
%This function is meant to be used on windows machines, although
%it may work on unix as well.
%
%This function will recursively copy all the files listed in the source
%directory to the DestDir.  Before copying, this function looks at the
%modification dates of files that are named the same to see if they need
%copying.  This makes it useful as an automatic backup tool, as it can be
%run say once a week and it will only copy the files that have been
%modified.  I have been using this program for awhile now and it works
%well, but no doubt that statement will be of little comfort if this script
%hoses your harddrive.  Since this deals with file transfers there is always
%a possibility that something can go wrong and overwrite files that the user
%did not want overwrote.  I strongly encourage anyone who uses this function
%to try it on a few test case directories to convince yourself that it works.
%Also, this script won't touch the source files other than for copying except
%in cases where the source files are hidden.  In that case the function will
%remove the system and hidden file attributes then attempt to copy, and if
%succesful will reset the hidden and system attributes when done copying.
%If this type of behavior is undesireable set the ChangeAttribIfHidden flag
%to 0 below.
%Note:  This function uses lots of for loops because it uses cell arrays,
%       which can make execution slow for lots of files.  On the order of a
%       few minutes to determine which files should be copied.  Also, I'm 
%       using a brute force approach to determine which files should be 
%       copied, which I'm sure is pretty inefficient.  So if you care to do
%       so, feel free to modify this function.
%
%
% BackUpSourceDirs = {
%                     'C:\Installfiles'
%                     'C:\Eudora5'
%                     'C:\Documents and Settings\Matt Becker\My Documents'
%                     'C:\Misc'
%                    };
% DestDir = 'D:';
%
%
% modified by Ryan H 2019

tic;
ChangeAttribIfHidden = 1;
SkipMessage = 1;    %initalizes flag to display warning message

               
NumFiles = 0;             %Initialize to zero
DirFlag = 1;              %Initialize Sub-Directories Exist Flag
%Check to make sure Source directories exist
disp('Analyzing Source Directories....');
for i=1:length(BackUpSourceDirs)
    if(exist(BackUpSourceDirs{i}) == 0)
        error(['The source directory ' BackUpSourceDirs{i} ' does not exist']);
    end
end
%Now Loop Over Each Source Directory and compile a list of all the files
for i=1:length(BackUpSourceDirs)
    Source{i} = dir(BackUpSourceDirs{i});
    if(length(Source{i}) <= 2)   %Check to see if the directory is empty
        error(['Source Directory ' BackUpSourceDirs{i} ' is empty']);
    end
end
%Now we need to do a recursive search through each source directory and if
%the contents are also directories, we need to do add them to the source
%directory list and perform a recursive search through them
StartDir = 1;      %Index for which directory to do recursive search on
while(DirFlag)
    BreakFlag = 0;
    NumNewDirs = 0;          %Initializes nuumber of new directoires to 0
    %Now Loop Over Each Source Directory and compile a list of all the files
    for i=1:length(BackUpSourceDirs)
        Source{i} = dir(BackUpSourceDirs{i});
    end
    %Now check to see if any contents are a directory
    for i=StartDir:length(BackUpSourceDirs)
        if(~BreakFlag)   %Stops if a directory was added to the list
            for j=1:length(Source{i})
                if(Source{i}(j).isdir == 1 && ~strcmp(Source{i}(j).name,'.') && ~strcmp(Source{i}(j).name,'..')) %indicates a directory
                    %Add directory to list
                    BackUpSourceDirs{end + 1} = [BackUpSourceDirs{i} filesep Source{i}(j).name];
                    NumNewDirs = NumNewDirs + 1;
                    BreakFlag = 1;   %don't add any more source directories this loop
                    StartDir = i+1;  %increment starting directory index
                end
            end
        end
    end
    if(NumNewDirs == 0)   %no new directories to add so exit loop
        DirFlag = 0;
    end
end
%
%Now Loop Over Each Source Directory one final time to get files
for i=1:length(BackUpSourceDirs)
    Source{i} = dir(BackUpSourceDirs{i});
    %If source directory is empty is doesn't need copied so set to
    %empty matrix for easy checking later on
    if(length(Source{i}) <= 2) 
        disp(['Warning! Source Directory ' BackUpSourceDirs{i} ' is empty and will not be copied']);
        BackUpSourceDirs{i} = [];
    end
end
disp('Analyzing Destination Directories....');
%By default, all Destination Directories will be identical
%to source directories except for destination path
%Check for trailing \ on DestDir and remove if it has them
if(strcmp(DestDir(end),filesep))
    DestDir(end) = [];
end
BackUpDestDirs = BackUpSourceDirs;    %initializes BackUpDestDirs
for i=1:length(BackUpDestDirs)
    if(~isempty(BackUpDestDirs{i}))
        BackUpDestDirs{i} = [DestDir BackUpSourceDirs{i}];
        BackUpDestDirs{i}(length(DestDir)+1:length(DestDir)+2) = [];
    end
end
%Now make sure each destination directory exits and create if not
for i=1:length(BackUpDestDirs)
    if(~isempty(BackUpDestDirs{i}))
        if(exist(BackUpDestDirs{i}) ~= 7)   %7 represents a directory
            if(~SkipMessage)
                disp(['Warning the Destination directory ' BackUpDestDirs{i} ' does not exist']);
                disp('The script is paused.  If this is okay ENTER, otherwise press CNTRL-C');
                disp('This message will be displayed for every directory that needs to be created');
                letter = input('To skip this message in the future press s then ENTER:  ','s');
                if(strcmpi(letter,'s'))
                    SkipMessage = 1;
                end
            end
            mkdir(BackUpDestDirs{i});
        end
    end
end
%Now Loop Over Each Destination Directory and compile a list of all the files
for i=1:length(BackUpDestDirs)
    if(~isempty(BackUpDestDirs{i}))
        Dest{i} = dir(BackUpDestDirs{i});
    end
end
disp('Determining which files need to be copied....');
%Now loop over each source directory and for each file check to see if its
%contained in the destination directory.  If not we need to copy it.  If
%so check modification dates.  If the source is newer copy otherwise don't
for i=1:length(BackUpSourceDirs)    %for each source directory
    if(~isempty(BackUpDestDirs{i}))  %as long as the directory is not empty
        for j=1:length(Source{i})   %for each file
            CopyFlag{i}(j) = 1;   %Initialize Copy Flag Assume Copy unless a match is found
            for(k=1:length(Dest{i}))  %for each destiniation file
                if(strcmp(Source{i}(j).name,Dest{i}(k).name))   %If files match assume don't copy
                    CopyFlag{i}(j) = 0;
                    if(~strcmp(Source{i}(j).date,Dest{i}(k).date))   %If dates don't match assume copy
                        CopyFlag{i}(j) = 1;
                    end
                end
                if(Source{i}(j).isdir == 1) %indicates a directory
                    CopyFlag{i}(j) = 0;    %don't copy directories
                end
            end
        end
    end
end
%Now copy all files indicated by a copyflag value of 1
for i=1:length(BackUpSourceDirs)
    if(~isempty(BackUpDestDirs{i}))
        for j=1:length(CopyFlag{i})
            if(CopyFlag{i}(j) == 1)
                if(strcmp(Source{i}(j).name,'.') || strcmp(Source{i}(j).name,'..'))
                    error(['Attempt to copy source file ' BackUpSourceDirs{i} filesep Source{i}(j).name]);
                end
                disp(['Copying ' BackUpSourceDirs{i} filesep Source{i}(j).name ' to ' BackUpDestDirs{i} filesep Source{i}(j).name]);
                try
                    copyfile([BackUpSourceDirs{i} filesep Source{i}(j).name],[BackUpDestDirs{i} filesep Source{i}(j).name]);
                    NumFiles = NumFiles + 1;
                catch
                    if(ChangeAttribIfHidden)
                        try  %try changing the file attribute to unhidden
                            fileattrib([BackUpSourceDirs{i} filesep Source{i}(j).name],'-s -h');
                            copyfile([BackUpSourceDirs{i} filesep Source{i}(j).name],[BackUpDestDirs{i} filesep Source{i}(j).name]);
                            NumFiles = NumFiles + 1;
                            fileattrib([BackUpSourceDirs{i} filesep Source{i}(j).name],'+h +s');
                        catch
                            disp(['Warning! An Error Occured when Copying ' BackUpSourceDirs{i} filesep Source{i}(j).name]);
                            disp('This file could not be copied!!');
                        end
                    else
                        disp(['Warning! An Error Occured when Copying ' BackUpSourceDirs{i} filesep Source{i}(j).name]);
                        disp('This file could not be copied!!');
                    end
                end
            end
        end
    end
end
TotalTime = toc;
TotalMin = floor(toc/60);
TotalSec = TotalTime - TotalMin*60;
disp([num2str(NumFiles) ' file(s) copied in ' num2str(TotalMin) ' mintue(s) and ' num2str(floor(TotalSec)) ' second(s).']);