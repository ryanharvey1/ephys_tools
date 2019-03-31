function [TC,Occ]= TuningCurvesBound(S, varargin)

% [TC,Occ] = TuningCurvesBound(S, X1, b1, min1, max1, X2, b2, min2, max2, ...)
% 
% INPUTS:
%      S -- a cell array of class ts or a single object of class ts
%      Xi, bi -- pairs of class ctsd or tsd (X) & a number (b)
%         Xi = the data along which we are finding the tuning (i.e. head direction "HD"
%              or position along track "V")
%         bi = number of bins in which to separate that dimension 
%              if b is omitted for a dimension, it defaults to 64
%         mini = minimum environment boundary in this dimension
%         maxi = maximum environment boundary in this dimension
% 
% OUTPUTS:
%    TC -- a cell array (if S was a cell array) or a single (if S is one object)
%          n-dimensional matrix of number of spikes occuring at each bin
%    Occ -- Occupancy within each bin
%
%    NOTE: TC is not normalized by default.  To normalize TC, use TC/(Occ + epsilon)
%
% ADR 1998
% version L5.0
% status: PROMOTED
%
% v 4.1 2 nov 1998 now takes out all NaNs.
% v 4.2 17 nov 98 now correctly returns same order of dimensions
% v 5.0 8 dec 99 now uses timestamps from first dim input (X1) for occupancy samples
% ZN added boundaries of the environment as an input, so that bins based on
% environment boundaries, not min and max position visited

%--------------------
% Unpack inputs
% V = cell array of tsd or ctsd of each dimension
% nV = array of bins for each dimension
%--------------------
vi=1; vc=1;
while vi <= length(varargin)
    % Xi
    V{vc} = varargin{vi};
    % bi
    if vi+1 <= length(varargin) && isa(varargin{vi+1},'double')
        nV(vc) = varargin{vi+1};
        vi = vi+1;
        % mini
        if vi+1 <= length(varargin) && isa(varargin{vi+1},'double')
            mV(vc)=varargin{vi+1};
            vi=vi+1;
            % maxi
            if vi+1 <= length(varargin) && isa(varargin{vi+1},'double')
                MV(vc)=varargin{vi+1};
                vi=vi+1;
            else
                MV(vc) = max(Data(V{vc}));
            end
        % mini
        else
            mV(vc) = min(Data(V{vc}));
        end
    % bi
    else
        nV(vc) = 64;
    end
    vi=vi+1;
    vc=vc+1;
end

%--------------------
% Size/type checks

if ~isa(S, 'cell')
   S = {S};
end

for iC = 1:length(S)
   if ~isa(S{iC}, 'ts')
      error(['S{', num2str(iC), '} is not a ts object.'])
   end
end

%--------------------
% parameters
nCells = length(S);
nD = length(V);
epsilon = 1e-100;

%--------------------
% restrict data to be in range for which we have 
% sufficient data

Tmin = -inf; 
Tmax = inf;
for iV=1:nD
   Tmin = max(Tmin, StartTime(V{iV}));
   Tmax = min(Tmax, EndTime(V{iV}));
end

VR = cell(size(V));
for iV=1:nD
   VR{iV} = Restrict(V{iV}, Tmin, Tmax);
end

SR = cell(size(S));
for iS=1:nCells
   SR{iS} = Restrict(S{iS}, Tmin, Tmax);
end

% % get mins & maxes
% 
% mV = zeros(size(V)); 
% MV = zeros(size(V)); 
% for iV=1:nD
%    mV(iV) = min(Data(V{iV}));
%    MV(iV) = max(Data(V{iV}));
% end

% prepare output

M = cell(nCells,1);

%--------------------
% Calculate fields
%--------------------

% foreach cell in the SpikeList
for iC = 1:nCells
   
   spikeTimes = Data(SR{iC});
   
   if isempty(spikeTimes)
      TC{iC} = squeeze(zeros([nV 1]));	% need to add dummy dimension if 1D
   else  
      dataValue = zeros(nD,length(spikeTimes));   % nD x |S{iC}| array for ndhist
      for vd=1:nD
         dataValue(vd,:) = Data(VR{vd}, spikeTimes)';  % VR is restricted V
      end
      [nanRows, nanCols] = find(isnan(dataValue));
      dataValue(:,nanCols) = [];
      TC{iC} = ndhist(dataValue, nV', mV', MV');
      % need to permute TC to return dimensions in correct order
      if (length(nD) > 1)
	TC{iC} = permute(TC{iC}, nD:-1:1);
      end
   end    
end

if length(TC) == 1
   TC = TC{1};
end

%--------------------
% Calculate Occupancy
%--------------------

timestamps = Range(VR{1}, 'ts');            % timing of first inputs
dataValue = zeros(nD, length(timestamps));
for vd = 1:nD
   dataValue(vd,:) = Data(VR{vd}, timestamps)';
end
[nanRows, nanCols] = find(isnan(dataValue));
dataValue(:,nanCols) = [];
if ~isempty(nanCols)
    disp(['NaN occupancy values exist at: ', num2str(nanCols)])
end

Occ = ndhist(dataValue, nV', mV', MV');
if (length(nD) > 1)
  Occ = permute(Occ, nD:-1:1);
end

Occ = Occ * DT(VR{1}) / 10000; % normalize for time

