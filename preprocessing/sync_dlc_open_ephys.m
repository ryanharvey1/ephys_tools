function sync_dlc_open_ephys(varargin)

p = inputParser;
p.addParameter('path',pwd);
p.parse(varargin{:});

path = p.Results.path;
[~,basename] = fileparts(path);

% Import Dat 

% Import timestamps 

% Save synced timestamps back to DLC 


end
