function varargout = WScatterFields(S, varargin)

% [SF1, SF2, ... SFn] = WScatterFields(S, X1, X2, ... Xn, parameters)
%
% INPUTS
%     S = cell array of cells
%     Xi = [c]tsd of data
% 
% OUTPUTS
%     SFi = cell array of [c]tsd of data when cells fired a spike
%
% In a WScatter field, the data is nSpikes by 2, column one contains the actual data
% (e.g. the position of the animal at the time of the spike) while column 2 contains 
% the weight, so that the WScatterField is correctly normalized.
% 
% PARAMETERS
%    NormalizeByTotalSpikes, default 1
%    NormalizeBySpeed, default 1
%
% Uses scatter fields

% ADR 1998
% version L1.0
% Status PROMOTED

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
if isa(varargin{iV}, 'char')
   varargin = varargin(iV:end);
else
   varargin = {};
end

% --------------------
% Parameters
NormalizeByTotalSpikes = 1;
NormalizeBySpeed = 1;
Extract_varargin;

nDim = length(X);
nCells = length(S);

%----------------------
for iD = 1:nDim
   
   % generate normal scatter field
   SF{iD} = ScatterFields(S, X{iD});
   if ~isa(SF{iD}, 'cell')
      SF{iD} = {SF{iD}};
   end
   for iC = 1:nCells
      TIME = Range(SF{iD}{iC}, 'ts');
      DATA = Data(SF{iD}{iC});
      WT = ones(size(DATA));
      
      % normalize by total spikes
      if NormalizeByTotalSpikes 
         WT = WT ./ length(DATA);
      end
      
      % normalize by running speed
      if NormalizeBySpeed
         speed = Speed(X{iD});
         WT = WT .* Data(speed, TIME);
      end
      
      SF{iD}{iC} = tsd(TIME, [DATA WT]);
      
   end
end

%------------------
for iV = 1: length(SF)
   if size(SF{iV}) == [1,1]
      varargout{iV} = SF{iV}{1};
   else
      varargout{iV} = SF{iV};
   end
end
