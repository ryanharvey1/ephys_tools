function StatusWarning(status, funcname)

% StatusWarning  Prints a status warning
%
% StatusWarning(status, funcname)
%
% INPUTS:
%       status - status message (string)
%       funcname - name of function to print message about (string)
% OUTPUTS:
%       (none)
%
% ADR 1998, version L4.0, last modified '98 by ADR

% status: PROMOTED


disp(['WARNING: "', funcname,'" is still under status --', status, '-- use at your own risk.']);
