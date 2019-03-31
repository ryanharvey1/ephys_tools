% MoveRawData
rats={'RH13','RH14','RH11','RH16','RH17','RH19','RH21','RH23','LS17','LS19','LS21','LS23','LE2813','LE2821','LE2823'};
for irats=1:length(rats)
    parent=strcat('F:\Users\reharvey\Place_Cell_Data\PAE_Rat\',rats(irats));
    disp(['CYCLING THROUGH RAT:',char(rats(irats))])
    parent=char(parent);
    structdir=dir(parent);
    for I=1:length(structdir)
        if structdir(I).isdir && structdir(I).name(1) ~= '.'
            if sum(ismember(structdir(I).name,'p'))==0
                if ~exist(['D:\Place_Cell_Data\RawPAE_PlaceCell\',rats{irats}],'dir')
                    mkdir(['D:\Place_Cell_Data\RawPAE_PlaceCell\',rats{irats}])
                end
                copyfile([parent filesep structdir(I).name],['D:\Place_Cell_Data\RawPAE_PlaceCell\',rats{irats},filesep structdir(I).name]) 
            end
            keep('I','K','structdir','parent','CurrentMat','rats','irats');
        end
    end
    disp(['DONE WITH:',char(rats(irats))])
end