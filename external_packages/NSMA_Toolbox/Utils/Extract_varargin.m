% Extract_varargin  (NOT A FUNCTION) Creates variables with names and values provided in varargin list of function
%
% Extract_varargin
%
% INPUTS:
%       (none)
% OUTPUTS:
%       (none)
%
% Allows use of varargin for parameter passing
% Expects varargin to consist of sequences of 'variable', value
% Sets variable to value for each pair.
% Changes the current workspace!
%
% ADR 1998, version L4.0, last modified '98 by ADR

% status: PROMOTED


for iV = 1:2:length(varargin)
    if ~ischar(varargin{iV})
        error('Extract_varargin:InputError','Variable input argument list must be a sequence of pairs (VariableNameString, Value)');
    end
    eval([varargin{iV}, ' = ', 'varargin{iV+1};']);
end