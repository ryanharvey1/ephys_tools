% CylinderDisplacement
% This script will cycle through a parent directory and pull out info from the .mat files
% input:
%       Laura=0; If you want Laura's data
%       Ryan=1; If you want Ryan's data
% output:
%       SAVES EXCEL SHEETS FOR EACH RAT WITH ALL VARIABLES
%
% Ryan Harvey 9/16
% CylinderDisplacement
tic
clear ; clc ; close all
addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analyses\spikeCode'));
addpath(genpath('d:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analyses\spikeCode'));

rats={'LS17','LS19','LS21','LS23','LE2813','LE2821'};
datafolder='F:\Users\reharvey\Place_Cell_Data\PAE_Rat\';

for irats=1:length(rats)
    parent=strcat(datafolder,rats(irats));
    disp(['CYCLING THROUGH RAT:',char(rats(irats))])
    parent=char(parent);
    structdir=dir(parent);
    K=1;
    for I=1:length(structdir) % 1 TO # OF FILES IN DIR
        if structdir(I).isdir && structdir(I).name(1) ~= '.' % IF DIR NOT '.'
            if exist(([parent filesep structdir(I).name filesep 'TT']),'file'); % LOCATE TT FOLDER
                cd([parent filesep structdir(I).name filesep 'TT']); % CD TO TT FOLDER
                CurentFolder=pwd;
                CurrentMat=dir('*.mat'); % LOCATE .MAT FILES IN TT FOLDER
                for J=1:length(CurrentMat) % 1 TO # OF .MAT FILES IN TT FOLDER
                    CurrentMatworking=CurrentMat(J).name;
                    if sum(ismember(CurrentMatworking,'spikeData'))==11 && sum(ismember(CurrentMatworking,'pathProperties'))~=16
                        sessionnum=ismember(CurrentMatworking,'S2');
                        if sum(sessionnum(1,end-15:end-14))~=2; continue; end
                        disp(['READING: ',CurentFolder filesep CurrentMatworking])
                        warning('off','MATLAB:load:variableNotFound');
                        split=strsplit(CurrentMatworking,'_');split{3}='S3';
                        try
                        CurrentMatworking2=strcat(split{1},'_',split{2},'__',split{3},'_',split{4});
                        catch 
                            continue
                        end
%                         load(CurrentMatworking,'nSpikes','OverallFR','FieldWidth','InformationContent');
%                         if nSpikes<50 || OverallFR>10 || FieldWidth<2 || InformationContent<0.3818; continue;end
                        load(CurrentMatworking,'SmoothRateMap');if ~exist('SmoothRateMap','var'); continue;end;m1=SmoothRateMap; 
                        if exist(CurrentMatworking2,'file')
                            load(CurrentMatworking2,'SmoothRateMap');if ~exist('SmoothRateMap','var'); continue;end;m2=SmoothRateMap; 
                        else
                            continue
                        end
%                         [D,M]=Displacement(SmoothRateMap1,SmoothRateMap2);
                            [d,c]=Displacement2(m1,m2);
%                             if d==270
%                                 test=1
%                             end
                        try
                            AllSpikeData(K,:)=[[CurentFolder filesep CurrentMatworking],num2cell(d),num2cell(c)];
                        catch
                            test=1
                        end
                        K=K+1; % CYCLE THROUGH OUTPUT FILE SO YOU DON'T OVERWRITE
                    end
                    if ~exist('AllSpikeData','var') % IF THE .MAT FILE WAS INCORRECT
                        keep('I','J','K','structdir','parent','datafolder','CurrentMat','CurrentFolder','rats','irats');
                    else
                        keep('I','J','K','AllSpikeData','structdir','parent','datafolder','CurrentMat','CurentFolder','rats','irats'); % KEEP ALLSPIKEDATA & STRUCTDIR IF .MAT WAS CORRECT
                    end
                end
            end
            if exist('AllSpikeData','var')
                keep('I','K','AllSpikeData','structdir','parent','datafolder','CurrentMat','rats','irats'); % DON'T KEEP VAR J BECAUSE YOU JUST EXITED THE J LOOP
            else
                keep('I','K','structdir','parent','datafolder','CurrentMat','rats','irats');
            end
        end
    end
if exist('AllSpikeData','var')
    AllSpikeData = cell2table(AllSpikeData,'VariableNames',{'CurentDir','DegreeDisplacement','DisplacementCorrelation'}); % CREATE TABLE
    [filepath, filename] = fileparts(parent);
    writetable(AllSpikeData,[filepath filesep '_CylinderDisplacement' filename '.xlsx']);
end

disp(['DONE WITH:',char(rats(irats))])
keep('rats','irats','parent','datafolder','filepath')

end

% ---------------------------------------------------------------------------------BRING BACK DATA TO ANALYZE-------------------------------------------------------------------------------------

% LS17 = xlsread([filepath,filesep,'_CylinderDisplacementLS17.xlsx']);
% LS19 = xlsread([filepath,filesep,'_CylinderDisplacementLS19.xlsx']);
% LS21 = xlsread([filepath,filesep,'_CylinderDisplacementLS21.xlsx']);
% LS23 = xlsread([filepath,filesep,'_CylinderDisplacementLS23.xlsx']);
% 
% control=[LS21;LS23]; 
% PAE=[LS17;LS19]; 
% numofcells=length(control)+length(PAE);
% clear LS17 LS19 LS21 LS23
% 
% [ AllStats ] = ScatterBox(control,PAE,{'Control' 'PAE'},{'DegreeDisplacement','DisplacementCorrelation'},2);
% 
% disp('done')
% toc




