function run_klusta(folder)
% Run klusta ephys pipeline for a given session folder across all shanks. 

% Find folders with klusta data 
[~,basename] = fileparts(folder);
d = dir(folder);
dfolders = d([d(:).isdir]);
klusta_folders = dfolders(contains({dfolders.name},'klusta'),:);

% iterate through and run klusta on each 
parfor shank = 1:length(klusta_folders)
    disp(['Running klusta on shank: ',num2str(shank)])
    path = [klusta_folders(shank).folder,filesep,klusta_folders(shank).name];
    command1 = ['activate klusta & cd ' path ' & klusta ' basename '_sh',num2str(shank),'.prm'];
    system(command1)
    
end
end