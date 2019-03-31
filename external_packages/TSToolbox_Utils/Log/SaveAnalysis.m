function SaveAnalysis(processdir,fname,variables,varNames,varargin)

% USAGE
%     SaveAnalysis(processdir,fname,variables,varNames,varargin)
%     
% this function saves analysis results in the 'Analysis' folder
% 
% INPUT:
%     processdir: the directory containing the Analysis folder. Typically, call it with 'pwd';
%     fname: filename of the mat files containing the analysis
%     variables: cell array of variable values
%     varNames: name of variables
%     
%     optionnal:
%      saveAnalysis(processdir,fname,variables,varNames,info)
%     info is cell array of char providing some essential info for each
%     variable (e.g. explaining the analysis or the format of the variable)
%     
%     Adrien Peyrache, 2012

if length(variables) ~= length(varNames)
    error('Variable values and variable names cell arrays must of the same size')
end

if ~isempty(varargin)
    info = varargin{1};
    if ~isa(info,'cell')
        error('Info must be a cell array')
    end
else
    info = cell(length(variables),1);
    for ii=1:length(info)
        info{ii} = 'No info';
    end
end

if length(info) ~= length(varNames)    
   error('Variable values and variable names cell arrays must of the same size')
end

analDir = [processdir filesep 'Analysis'];
if ~exist(analDir,'dir')
    mkdir(analDir)
end

append=1;

analFile = [analDir filesep fname];
if ~strcmp(analFile(end-3:end),'.mat')
    analFile = [analFile '.mat'];
end
if ~exist(analFile,'file')
    append=0;
end
    
for ii=1:length(variables)
    
    eval([varNames{ii} ' = variables{ii};'])
    eval([varNames{ii} '_Info = info{ii};'])

    if ii==1 && ~append
        save(analFile,varNames{ii})
        append = 1;
    else
        save(analFile,varNames{ii},'-append');
    end
   
    save(analFile,[varNames{ii} '_Info'],'-append');
    
end