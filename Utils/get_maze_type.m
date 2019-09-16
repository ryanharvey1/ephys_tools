function mazetype=get_maze_type(data)
% get_maze_type: gets maze type from rat's metadata file
%
% ryan harvey 2019

path_parts=strsplit(data.session_path,filesep);
path_parts(end-2:end)=[];
metadata_file=fullfile(strjoin(path_parts,filesep),'AnimalMetadata',[data.rat,'_metadata.mat']);
if exist(metadata_file,'file')
    load(metadata_file,'AnimalMetadata');
    sessions=fieldnames(AnimalMetadata.RecordingLogs);
    mtypes=strsplit(AnimalMetadata.RecordingLogs.(sessions{ismember(sessions,data.sessionID)}).MazeTypes,',');
    
    if isempty(mtypes{1})
        return % will result in an error
    end
    
    for mt=1:size(data.events,2)
        mazetype{mt}=mtypes{mt};
    end
    
    mazetype=strtrim(mazetype);
end