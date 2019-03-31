function varargout = ScatterFields(S, varargin)

% [SF1, SF2, ... SFn] = ScatterFields(S, X1, X2, ... Xn, parameters)
%
% INPUTS
%     S = cell array of cells
%     Xi = [c]tsd of data
% 
% OUTPUTS
%     SFi = cell array of [c]tsd of data when cells fired a spike
%
% Warning: Removed spike pair functionality.  Use SelectSpikePairs instead.
% Warning: Removed TrialMask flag.  Use mask directly.
% Warning: No longer skips cells, fills with blank tsd instead
%
% NOTE: To normalize scatter fields by total spikes, divide by len(Data(SF))
%       To normalize each pt by occupancy, normalize by running speed at time of each spike

% ADR 1998
% version L4.4
% Status PROMOTED

% v4.4 17 nov 98 no longer skips cells.

% --------------------
% Prepare IO
SF = {}; X = {};

% --------------------
% separate out varargin
for iV = 1:length(varargin)
   switch class(varargin{iV})
   case {'ctsd','tsd'}
      X{iV} = varargin{iV};
   case 'char'
      break;
   otherwise
      error('Unknown input class.');
   end
end
varargin = {};

% --------------------
% Parameters
Extract_varargin;

nDim = length(X);
nCells = length(S);

% --------------------
% Generate fields
for iD = 1:nDim
   for iC = 1:nCells
      varargout{iD}{iC} = Restrict(X{iD}, []);  % if is empty fill with blank space
      if ~isempty(S{iC}) 
         XS = Restrict(X{iD}, Data(S{iC}));
         r = Range(XS,'ts');
         d = Data(XS);
         f = find(~isnan(d));
         if ~isempty(f)
            varargout{iD}{iC} = tsd(r(f),d(f));
         end % if iC fires during trials
      end % if iC contains data
   end % for all cells
end % for all inputs

%------------------
for iV = 1: length(varargout)
   if size(varargout{iV}) == [1,1]
      varargout{iV} = varargout{iV}{1};
   end
end
