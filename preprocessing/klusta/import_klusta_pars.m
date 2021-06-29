function klustapars = import_klusta_pars(varargin)
%IMPORT_KLUSTA_PARS  Read parameters from defaults text file as a matrix.
%
%  klustapars = import_klusta_pars;
%
% By: Max Murphy  v1.0  01/03/2018  Original version (R2017b)

%% DEFAULTS
DELIM = {''};
FNAME = 'klusta_pars.prm';
FORMATSPEC = '%s%[^\n\r]';
STARTROW = 0;

%% PARSE VARARGIN
for iV = 1:2:numel(varargin)
   eval([upper(varargin{iV}) '=varargin{iV+1};']);
end

%% Open the text file.
fileID = fopen(FNAME,'r');

%% Read columns of data according to the format.
textscan(fileID, '%[^\n\r]', STARTROW, ...
   'WhiteSpace', '', ...
   'ReturnOnError', false);
dataArray = textscan(fileID, FORMATSPEC, inf,...
   'Delimiter', DELIM, ...
   'TextType', 'string', ...
   'ReturnOnError', false, ...
   'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Create output variable
klustapars = [dataArray{1:end-1}];


end
