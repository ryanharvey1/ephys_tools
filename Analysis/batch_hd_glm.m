function batch_hd_glm(dataset,save_path,sessions)
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

addpath(genpath('D:\Users\BClarkLab\ephys_tools\external_packages\ln-model-of-mec-neurons_laura'))
cd(dataset)

% Get cell id from compiled table
sessions = readtable(sessions);
hd_met = zeros(length(sessions.SessionID),1);
for i = 1 : length(sessions.SessionID)
    
    % check to see if the cell was already run. If so, continue
    if exist([save_path,'HD_results_',extractBefore(sessions.SessionID{i},'.mat'),'_',extractBefore(sessions.tetrode{i},'.mat'),'_',num2str(sessions.cell(i)),'.mat'],'file')
        continue
    end
    
    %Load Data session cell
    data = load([dataset,filesep,sessions.SessionID{i}]);
    sess = 1; % Lets just look at baseline 1 for now
    tet = sscanf(sessions.tetrode{i},'TT%d.mat');
    [celln] = find_cells(data,tet,sessions.cell(i));
    
    % Run GLM 
    disp('Running GLM')
    glm_res = glm_HD(data,sess,celln);
    
    % Run Shuffling procedure on cell to verify GLM results 
    disp('Running Shuffling Procedure')
    [shuff_pass,p_val,z]= shuff({sessions.SessionID{i},sessions.tetrode{i},sessions.cell(i)},'feature',{'dic'},'nshuffle',1000);
    
    glm_res.shuffling.shuff_pass = shuff_pass;
    glm_res.shuffling.p_val = p_val;
    glm_res.shuffling.z_score = z; 
    
    % Add binary indicator variable for ease
    if glm_res.best_model == 1 && p_val(1) < 0.05 
        glm_res.hd_met = 1;
        hd_met(i,1) = 1;
    else 
        glm_res.hd_met = 0;
        hd_met(i,1) = 0;
    end
    
    % Save file
    disp('saving results')
    save([save_path,'HD_results_',extractBefore(sessions.SessionID{i},'.mat'),'_',extractBefore(sessions.tetrode{i},'.mat'),'_',num2str(sessions.cell(i)),'.mat'],'glm_res')
end

% Update 
disp([num2str(sum(hd_met)),' HD cells found out of ', num2str(length(sessions.SessionID)), ' total cells.'])
% Save HD cell list
writetable(sessions(logical(hd_met),:),[save_path,'hd_cell_list.csv'])

end

