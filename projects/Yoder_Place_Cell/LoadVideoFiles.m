% LOAD VIDEO FILES
% MeanRateMaps
clear ; close all ; clc
addpath(genpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/NSMA_Toolbox'));
addpath(genpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis/Yoder_Place_Cell'));
addpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis');
addpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/CircStat2012a')
addpath('/Users/ryanharvey/GoogleDrive/MatlabDir/CStr2String')
% FigureLocation='/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/Figures_For_Paper';
% paths=importdata('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/GW- Field_StatsPlace Cells_Tilted_Mice_1_1.xlsx');
load('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/GW_Field_StatsPlaceCells_Tilted_Mice.mat')
% load('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/NewData.mat')
path_group=data.textdata.Field_StatsPlaceCells0x2DTilted(:,1:2);
path_group(1,:)=[];
group=path_group(:,2);
control_paths=path_group(strcmp(group,'Control'),1);
tilted_paths=path_group(strcmp(group,'Tilted'),1);


% just first session
Sess1con=control_paths(1:5:length(control_paths));
Sess1tilt=tilted_paths(1:5:length(tilted_paths));
% MAKE PATHS UNIQUE BEFORE RUNNING

% for ii=1:length(tilted_paths)
%     split=strsplit(tilted_paths{ii},'ASCII');
%     Sess1conW(ii,:)=split(1);
% end
% unique(Sess1conW(:,1));

Controlresults=[];
for i=1:length(Sess1con)
    try
        % Import video file
        split=strsplit(Sess1con{i},'ASCII');
        disp(split(1))
        fileID=fopen(char(strcat(split(1),'VT1_DupRem_PrcFields_lights.txt')),'r');
        Intro = textscan(fileID,'%s%s%s');
        fclose(fileID);
        
        X = [Intro{2}];
        Y = [Intro{3}];
        
        X=X(~cellfun('isempty',X));
        Y=Y(~cellfun('isempty',Y));
        
        VIDEO=regexprep([X(4:end,1),Y(3:end,1)],',','','emptymatch');
        N = sscanf(CStr2String(VIDEO, '*'), '%f*');
        B = reshape(N,[],2);
        

        % interp and smooth
%         video(isnan(video))=0;
%         [x,y] = InterpolateTrackerZeros_v2(video(:,2),video(:,3));
        smoothVideo=[SmoothVector(B(:,1),3),SmoothVector(B(:,2),3)];
        
        % calculate velo
        [vel_cmPerSec,vel_abs,pixelDist] = InstaVel(smoothVideo,'no',61,60);
        vel_cmPerSec(vel_cmPerSec>70)=[]; vel_abs(vel_abs>70/pixelDist/60)=[];
        
        ExtractedAngle=XYangle(smoothVideo(:,1),smoothVideo(:,2));
        
        % EXTRACT BASIC MOVEMENT DATA
        AverageAnglePerSec=rad2deg(circ_mean((abs(diff(deg2rad(ExtractedAngle))))*60));
        OverallDistTraveled=sum(vel_abs*pixelDist);
        MeanVelocity=mean(vel_cmPerSec);
        STDVelocity=std(vel_cmPerSec);
        medVelocity=median(vel_cmPerSec);
        
        Controlresults=[Controlresults;AverageAnglePerSec,OverallDistTraveled,MeanVelocity,STDVelocity,medVelocity];
    catch
        disp('ERROR SOMEWHERE...DARN IT''S OKAY')
    end
end

Tiltedresults=[];
for i=1:length(Sess1tilt)
    try
        % Import video file
        split=strsplit(Sess1tilt{i},'ASCII');
        disp(split(1))
        fileID=fopen(char(strcat(split(1),'VT1_DupRem_PrcFields_lights.txt')),'r');
        Intro = textscan(fileID,'%s%s%s');
        fclose(fileID);
        
        X = [Intro{2}];
        Y = [Intro{3}];
        
        X=X(~cellfun('isempty',X));
        Y=Y(~cellfun('isempty',Y));
        
        VIDEO=regexprep([X(4:end,1),Y(3:end,1)],',','','emptymatch');
        N = sscanf(CStr2String(VIDEO, '*'), '%f*');
        B = reshape(N,[],2);
        

        % interp and smooth
%         video(isnan(video))=0;
%         [x,y] = InterpolateTrackerZeros_v2(video(:,2),video(:,3));
        smoothVideo=[SmoothVector(B(:,1),3),SmoothVector(B(:,2),3)];
        
        % calculate velo
        [vel_cmPerSec,vel_abs,pixelDist] = InstaVel(smoothVideo,'no',61,60);
        vel_cmPerSec(vel_cmPerSec>70)=[]; vel_abs(vel_abs>70/pixelDist/60)=[];
        
        ExtractedAngle=XYangle(smoothVideo(:,1),smoothVideo(:,2));
        
        % EXTRACT BASIC MOVEMENT DATA
        AverageAnglePerSec=rad2deg(circ_mean((abs(diff(deg2rad(ExtractedAngle))))*60));
        OverallDistTraveled=sum(vel_abs*pixelDist);
        MeanVelocity=mean(vel_cmPerSec);
        STDVelocity=std(vel_cmPerSec);
        medVelocity=median(vel_cmPerSec);
        
        Tiltedresults=[Tiltedresults;AverageAnglePerSec,OverallDistTraveled,MeanVelocity,STDVelocity,medVelocity];
    catch
        disp('ERROR SOMEWHERE...DARN IT''S OKAY')
    end
end


% VarNames=regexprep(VarNames,'R','','emptymatch');
%
% fid = fopen(filename,'r');
% MyTextFile = textscan(fid,'%s','delimiter',',');
% fclose(fid);
% MyTextFile = [MyTextFile{:}];
%
% MyTextFile=MyTextFile(9:end,1);
% MyTextFile(strcmp(MyTextFile,'R'))=[];
%
% eval(MyTextFile{2});
% eval(MyTextFile{3});
% % extract the part containg the matrix to read
% MatrixLines = MyTextFile(9:end);
% % get rid of Chirp1, Chirp2, ...
% MatrixLines = regexprep(MatrixLines,'Chirp[^%d]','');
% % Convet the char into a numerical matrix.
% MyMatrix = str2num(cell2mat(MatrixLines));
% %
%
% fid   = fopen(filename);
% line2 = textscan(fid, '%f%f%f%f%f%f%f%f%f\r\n %*[^\n]','HeaderLines',6);
% fclose(fid);



% MyTextFile=mat(MyTextFile);
%
% MyTextFile{MyTextFile>1000}=[];



