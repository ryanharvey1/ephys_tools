function handleAnimalMetaData(basepath,rat)
% handleAnimalMetaData: creates and updates animal meta data
% This function can handle most, if not all, of your animal and session meta
% data. Upon calling, this function will initally prompt you for specific information
% about your animal. It will then save this info in the specified folder.
% You can then call the function again and it will ask you about your
% recording sessions. All info is again saved to the same specified folder.
%
% Input:
%       basepath: path where your animal meta data will live (should be close to your data)
%       rat: animal ID
%       
%       you can also call this function with no inputs and it will prompt
%       you for what it needs
%
% Output:
%       no specific outputs, but this function will save your metadata .mat
%       file
%
%
% example:
%
% >> handleAnimalMetaData('D:\Projects\PAE_PlaceCell\PAE_HPC_DATA\AnimalMetadata','LEM3116')
% Loading Metadata file for LEM3116
% Path to Session D:\Projects\PAE_PlaceCell\PAE_HPC_DATA\LEM3116\2018-07-13_11-12-38
% MazeTypes  linear track, cylindar,cylindar
% Number of Turns  1.5
% View current coordinates over brain atlas? ('y','n') y
% RecordingArea  ca1
% Notes  very good movement!
%  ---------------------------------------------------------
% More Sessions to input? ("y" or "n")  n
% Saving Infomation...
% done
%
% Ryan E Harvey

com=which('handleAnimalMetaData');
com=strsplit(com,filesep);
basedir=[com{1},filesep,'Users',filesep,com{3},filesep,'GoogleDrive',filesep,'MatlabDir'];
addpath([basedir,filesep,'BClarkToolbox',filesep,'Analysis',filesep,'Visualize'])

if ~exist('basepath','var')
    basepath=uigetdir(cd,'Locate Metadata folder');
    rat=input('Animal ID  ','s');
end

% Create stucture if metadata file does not exist
if exist([basepath,filesep,rat,'_metadata.mat'],'file')~=2
    % basic rat info
    AnimalMetadata.AnimalName = rat;
    AnimalMetadata.Animal.Species = input('Species  ','s');
    AnimalMetadata.Animal.Strain = input('Strain  ','s');
    AnimalMetadata.Animal.GeneticLine = input('GeneticLine  ','s');
    AnimalMetadata.Animal.Sex = input('Sex  ','s');
    AnimalMetadata.Animal.DateOfBirth = input('DateOfBirth YYYYMMDD  ');%YYYYMMDD format
    AnimalMetadata.Animal.WeightGramsAtSurgery = input('WeightGramsAtSurgery  ');%grams
    AnimalMetadata.Animal.ExperimentalTreatment = input('ExperimentalTreatment (PAE, Sac, etc.)  ','s');%grams
    
    % surgery info
    AnimalMetadata.Surgery.Date = input('Surgery Date YYYYMMDD  ');
    AnimalMetadata.Surgery.Anesthesia = 'Isoflurane';
    AnimalMetadata.Surgery.Analgesic = 'Buprenex';
    AnimalMetadata.Surgery.Complications = input('Surgery Complications  ','s');
    AnimalMetadata.Surgery.Notes = input('Surgery Notes  ','s');
    
    % manipulation info
    AnimalMetadata.LesionorInactivation.status = input('Did you cause a Lesion or Inactivation (''y'',''n'')  ','s');
    if contains(AnimalMetadata.LesionorInactivation.status,'y')
        AnimalMetadata.LesionorInactivation.TargetRegions=input('TargetRegion  ','s');
        AnimalMetadata.LesionorInactivation.type=input('type of lesion (muscimol,NMDA,etc)  ','s');
        AnimalMetadata.LesionorInactivation.Coordinates.Anteroposterior=input('Anteroposterior  ');
        AnimalMetadata.LesionorInactivation.Coordinates.Mediolateral=input('Mediolateral  ');
        AnimalMetadata.LesionorInactivation.Coordinates.DorsalVentral=input('DorsalVentral  ');
    end
    
    % ephys info
    AnimalMetadata.ExtracellEphys.Probes.mmPerScrewTurn = input('mmPerScrewTurn  ');
    AnimalMetadata.ExtracellEphys.Probes.NumberOfProbes = input('NumberOfProbes  ');
    AnimalMetadata.ExtracellEphys.Probes.TargetRegions = input('TargetRegions  ','s');
    AnimalMetadata.ExtracellEphys.Probes.TargetHemisphere = input('TargetHemisphere  ','s');
    for i=1:AnimalMetadata.ExtracellEphys.Probes.NumberOfProbes
        AnimalMetadata.ExtracellEphys.Probes.ImplantCoordinates.Anteroposterior(i,1) = input('Anteroposterior  ');%one for each probe
        AnimalMetadata.ExtracellEphys.Probes.ImplantCoordinates.Mediolateral(i,1) = input('Mediolateral  ');
        AnimalMetadata.ExtracellEphys.Probes.ImplantCoordinates.DorsalVentral(i,1) = input('DorsalVentral  ');
    end
    AnimalMetadata.ExtracellEphys.Probes.ImplantAngle.Anteroposterior = input('ImplantAngle.Anteroposterior  ');%degrees of top anterior as sitting behind animal
    AnimalMetadata.ExtracellEphys.Probes.ImplantAngle.Mediolateral = input('ImplantAngle.Mediolateral  ');%degrees clockwise as sitting behind animal
    AnimalMetadata.ExtracellEphys.Channels.ImpedanceFilenames = input('ImpedanceFilenames  ','s');%Filenames in basepath folder, or leave as {} if none
    
    disp('Saving Infomation to .mat...')
    save([basepath,filesep,rat,'_metadata.mat'],'AnimalMetadata')
    disp('done')
    return
end

% update recording logs
disp(['Loading Metadata file for ',rat])
load([basepath,filesep,rat,'_metadata.mat'],'AnimalMetadata')
moresessions='y';
while contains(moresessions,'y')
    if isfield(AnimalMetadata,'RecordingLogs')==0
        sessiondate=input('Path to Session ','s');
        ratID=strsplit(sessiondate,filesep);
        sessiondate=['S',strjoin(regexp(ratID{end},'\d*','Match'),'')];
        AnimalMetadata.RecordingLogs.(sessiondate).MazeTypes=input('MazeTypes  ','s');
        for i=1:length(AnimalMetadata.ExtracellEphys.Probes.ImplantCoordinates.Anteroposterior)
            AnimalMetadata.RecordingLogs.(sessiondate).DorsalVentral=input('Number of Turns  ')...
                *AnimalMetadata.ExtracellEphys.Probes.mmPerScrewTurn+...
                AnimalMetadata.ExtracellEphys.Probes.ImplantCoordinates.DorsalVentral;
        end
        if contains(input('View current coordinates over brain atlas? (''y'',''n'') ','s'),'y')
            brainAtlas([AnimalMetadata.ExtracellEphys.Probes.ImplantCoordinates.Mediolateral,...
                AnimalMetadata.ExtracellEphys.Probes.ImplantCoordinates.Anteroposterior,...
                AnimalMetadata.RecordingLogs.(sessiondate).DorsalVentral])
        end
        AnimalMetadata.RecordingLogs.(sessiondate).RecordingArea=input('RecordingArea  ','s');
        AnimalMetadata.RecordingLogs.(sessiondate).Notes=input('Notes  ','s');
    else
        disp('Sessions')
        % reorder RecordingLogs fields by ASCII order
        AnimalMetadata.RecordingLogs=orderfields(AnimalMetadata.RecordingLogs);
        disp(AnimalMetadata.RecordingLogs)
        disp('---------------------------------------------------------')
        sessiondate=input('Path to Session ','s');
        ratID=strsplit(sessiondate,filesep);
        sessiondate=['S',strjoin(regexp(ratID{end},'\d*','Match'),'')];
        AnimalMetadata.RecordingLogs.(sessiondate).MazeTypes=input('MazeTypes  ','s');
        sessions=fieldnames(AnimalMetadata.RecordingLogs);
        for i=1:length(AnimalMetadata.ExtracellEphys.Probes.ImplantCoordinates.Anteroposterior)
           try
               AnimalMetadata.RecordingLogs.(sessiondate).DorsalVentral(i,1)=input('Number of Turns  ')...
                *AnimalMetadata.ExtracellEphys.Probes.mmPerScrewTurn+...
                AnimalMetadata.RecordingLogs.(sessions{end-1}).DorsalVentral(i);
           catch
               AnimalMetadata.RecordingLogs.(sessiondate).DorsalVentral(i,1)=0;
           end
        end
        if contains(input('View current coordinates over brain atlas? (''y'',''n'') ','s'),'y')
            brainAtlas([AnimalMetadata.ExtracellEphys.Probes.ImplantCoordinates.Mediolateral,...
                AnimalMetadata.ExtracellEphys.Probes.ImplantCoordinates.Anteroposterior,...
                AnimalMetadata.RecordingLogs.(sessiondate).DorsalVentral])
        end
        AnimalMetadata.RecordingLogs.(sessiondate).RecordingArea=input('RecordingArea  ','s');
        AnimalMetadata.RecordingLogs.(sessiondate).Notes=input('Notes  ','s');
    end
    disp('---------------------------------------------------------')
    moresessions=input('More Sessions to input? ("y" or "n")  ','s');
    
    % reorder RecordingLogs fields by ASCII order
    AnimalMetadata.RecordingLogs=orderfields(AnimalMetadata.RecordingLogs);
end

disp('Saving Infomation...')
save([basepath,filesep,rat,'_metadata.mat'],'AnimalMetadata')
disp('done')
end