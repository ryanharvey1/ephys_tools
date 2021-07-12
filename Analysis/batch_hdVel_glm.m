function batch_hdVel_glm(dataset,save_path,sessions)
% Runs through ephys_tools processed data files to perform Hardcastle LN GLM (HD) on cells
% that have >1 Hz peak rate and >100 spikes. Dependent on ephys_tools
% formatting and toolbox contents. Requires:
%   - external_packages > ln-model-of-mec-neurons_laura
%   - external_packages > CircStat2012a
%   - Analysis > compile_cell_measures.m
%   - Utils > shuff
%
% Input:
%   dataset: path to processed data folder
%   save_path: path to save GLM results
%   sessions: path where list of cells to run (.csv output of
%   compile_cell_measures) resides
%
% Saves GLM and Shuffling results as mat files in in save_path. Saves
% updated csv for all cells that met both HD for GLM and for mean vector
% length and directional information content shuffling.
%
% Laura B 10/2020

% Get cell id from compiled table
sessions = flipud(readtable(sessions));
for i = 1 : length(sessions.SessionID)
    
    %check to see if the cell was already run. If so, continue
    if exist([save_path,'HDVel_results_',extractBefore(sessions.SessionID{i},'.mat'),'_',extractBefore(sessions.tetrode{i},'.mat'),'_',num2str(sessions.cell(i)),'.mat'],'file')
        disp('GLM completed. Running next cell.')
        continue
    elseif exist([save_path,'HDVel_results_',extractBefore(sessions.SessionID{i},'.mat'),'_',extractBefore(sessions.tetrode{i},'.mat'),'_',num2str(sessions.cell(i)),'.mat'],'file')
        disp('GLM completed. Running next cell.')
        continue
    end
    
    %Load Data session cell
    data = load([dataset,filesep,sessions.SessionID{i}],'frames','spikesID','maze_size_cm','samplerate','Spikes','events');
    sess = 1; % Lets just look at baseline 1 for now
    tet = sscanf(sessions.tetrode{i},'TT%d.mat');
    [celln] = find_cells(data,tet,sessions.cell(i));
    % Run GLM
    disp('Running GLM')
    glm_res = glm_HDVel(data,sess,celln);
    % Save file
    disp('saving results')
    save_glm_res(save_path,i,sessions,glm_res)
    
end

end

function save_glm_res(save_path,idx,sessions,glm_res)
save([save_path,'HDVel_results_',extractBefore(sessions.SessionID{idx},'.mat'),'_',extractBefore(sessions.tetrode{idx},'.mat'),'_',num2str(sessions.cell(idx)),'.mat'],'glm_res');
end

