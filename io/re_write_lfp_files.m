% re-write all .lfp files to 1250hz

% find all lfp files
lfp_files = dir('F:\Projects\PAE_PlaceCell\data\**\*.lfp');

Fold = 32000;
Fnew = 1250;
WaitMessage = parfor_wait(length(lfp_files),'Waitbar',true);
% loop through each file
for i = 1:length(lfp_files)
    
    % load raw data
    lfpfile = locate_ncs(lfp_files(i).folder);
    % preallocate & get time stamps
    [ts] = Nlx2MatCSC(lfpfile{1},[1 0 0 0 0], 0, 1, [] );
    signal = zeros(length(lfpfile),(length(ts)*512)*Fnew/Fold);
    
    % resample
    [fold, fnew] = rat(Fold./Fnew);
    for ii=1:length(lfpfile)
        fprintf(' %d ', ii);
        
        [~, filename] = fileparts(lfpfile{ii});
        ch = str2double(extractAfter(filename,'CSC'));
        
        [Samples]= Nlx2MatCSC(lfpfile{ii}, [0 0 0 0 1], 0, 1);
        
        signal(ch,1:length(Samples(:))*fnew/fold) =...
            resample(Samples(:), fnew, fold);
    end
    
    % save .lfp
    disp('saving .lfp file');
    fidout = fopen(fullfile(lfp_files(i).folder,lfp_files(i).name), 'w');
    fwrite(fidout,signal,'int16');
    fclose(fidout);
    WaitMessage.Send;
    
end
WaitMessage.Destroy;

function lfpfile = locate_ncs(session_path)
channels=table2cell(struct2table(dir([session_path,filesep,'*.ncs'])));
lfpfile=strcat(channels(:,2),filesep,channels(:,1));
[~,idx]=sort(str2double(extractBetween(lfpfile,'CSC','.ncs')));
lfpfile = lfpfile(idx);
end