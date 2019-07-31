function [thetaout]=fixNLXangle(thetain,max_allowed_flips)
% Corrects errors in headdirection tracking
%
% Input: 
%       thetain: head angle in degrees
%
%       recursion: 1 or 0 if you would like to smooth out all
%           non-detects... if 0, the max allowed flips or the maximum number of
%           consecutive lost samples is defaults to 5
% Output: 
%       thetaout: corrected head angle -180 to 180

% Neuralynx head direction recordings often have missed samples, registered
% as either 0 or NaN. This script finds all missed samples and accidental
% flips of 180 degree (or any large flips > 162 degress) that persist for more than int
% 'max_allowed_flips' (5 samples) and interpolates the directional data
% in-between. Then we smooth the vector by convolving it with a gaussian of
% standard deviation of 3 bins, (or .1 sec at 30hz sampling freq).
%
% For all numerical samples (0 excluded) we continuize the data so
% that there are no more jumps at 360 and 0, and linearly interpolate all
% areas of bad data (either large angular velocities or 450s, NaNs, or
% zeros) that are 'max_allowed_flips' or less. Then smooth, and bring back
% to [-180,180).
%

bads = (thetain==0 | thetain==450);
thetain(bads) = NaN;

% circular derivative
thetain=pi*thetain./180; thetaout = thetain; % convert thetain to radians

delta_theta=diff(thetain(~isnan(thetain))); % take difference in discontinuous angles in radius

delta_theta=atan2(sin(delta_theta), cos(delta_theta)); % derivative of circular data

thetaout(~isnan(thetaout)) = [thetain(find(~isnan(thetain),1,'first'));...
    cumsum(delta_theta) + thetain(find(~isnan(thetain),1,'first'))]; % continuized data
    
% if recursion==1
%     [start,ends,~]=findgroups(isnan(thetaout));
%     thetaout=RemoveJumpsAndSmooth(thetaout,max(ends-start)+1,.5*pi);
% else
    thetaout=RemoveJumpsAndSmooth(thetaout,max_allowed_flips,.5*pi);
% end

% fix the rest of NaN values
% thetaout=fillmissing(thetaout,'pchip',2,'EndValues','nearest');

thetaout=atan2(sin(thetaout),cos(thetaout)); % bring back to radians bound by [-pi, pi)

thetaout=(thetaout*180)./pi;

if iscolumn(thetaout)
    thetaout=thetaout';
end
end

function x = RemoveJumpsAndSmooth(x, max_allowed_flips, jitter_threshold, smooth_std,filter_def)
%  Removes large jumps in data, smoothes the removed regions
%   ARGUMENTS
%       x                   a vector of something. NaNs may be used in place of
%                           known bad samples
%       max_allowed_flips   (samples) How many consecutive samples can be
%                           corrected? default = 5
%       jitter_threshold    (units of x) What intersample interval
%                           constitutes a jump?
%                           default = 15
%       smooth_std          (samples) STD of gaussian smoothing kernel.
%                           default = 2
%
% This function takes a vector, x, and looks for jumps greather than
% jitter_threshold that last for less than max_allowed_flips in samples. We
% then linearly interpolate the missed samples and smooth by convolution
% with a gaussion of standard deviation smooth_std bins.
%
% x = RemoveJumpsAndSmooth(x, max_flip_length, jitter_threshold, smooth_std);

%     import CMBHOME.Utils.*

x = x(:);

ts = 0:length(x)-1;

ts2 = ts;

if ~exist('jitter_threshold', 'var'), jitter_threshold = 15; end % pixels; one-sample change in distance that qualifies as bad vector

if ~exist('max_allowed_flips', 'var'), max_allowed_flips = 5; end % samples

if ~exist('smooth_std', 'var'), smooth_std = 2; end % samples

if ~exist('filter_def','var'), filter_def = -6:6;end

% there are two cases for which are the flips, and which is the real data. we will see how many
% samples fall into each case, and then pick the lesser. then we can
% cut out any jumps that are too long

flips = findOnsetsAndOffsets(isnan(x));
if ~isempty(flips)
    flips(:,2) = flips(:,2)+1;
    
    flips = [flips(:); find(abs(diff(x))>jitter_threshold)+1]; % all flip points
    
    flips = sort(unique(flips)); % indeces of NaNs or jumps
    
    flips = [flips(1:end-1), flips(2:end)];  % epochs formation
    
    flips(:,2) = flips(:,2)-1; % adjust for diff shift
    
    flips(flips(:,2)-flips(:,1)>max_allowed_flips-1,:) = []; % remove those greater than max_allowed_flips
    
    flips = mat2cell(flips, ones(size(flips,1),1),2); % convert to pairs corresponding to steps
    
    flips = cellfun(@(c) c(1):c(2), flips, 'unif', 0); % convert to indices
    
    x([flips{:}]) = []; % remove samples in ts and x
    ts([flips{:}]) = [];
    
    x = interp1(ts, x, ts2);
end
% gaussian_filter = pdf('Normal',-5*smooth_std:5*smooth_std,0,smooth_std);
if smooth_std>0
    x = ndnanfilter(x, normpdf(filter_def, 0, smooth_std)', [], 1, {}, {}, 1); % convolve with gaussian and ignore NaNs
end
% x = conv(x,gaussian_filter, 'same'); % smooth by convolution with a gaussian

end

function [OnOffs] = findOnsetsAndOffsets(boolVec)
% Returns list of aligned start and stops of chunks of 1's in a vector of
% 1's and 0's.
%
% INPUTS
%  boolVec - vector of 1's and 0's
%
% OUTPUTS
%  startEnds - Nx2 list of indices of the first and last 1's for the N
%              contiguous blocks of 1's.
%
% % function [OnOffs] = findOnsetsAndOffsets(boolVec)


boolVec = boolVec(:)';

starts = find(diff(boolVec)==1);
ends = find(diff(boolVec)==-1);

% if starts out going correct speed, add 1 to starts
if boolVec(1)
    starts = [0 starts];
end

% if finishes going correct speed, add final value to ends
if boolVec(end)
    ends = [ends length(boolVec)];
end

OnOffs = [starts(:)+1 ends(:)];
end
function [Y,W] = ndnanfilter(X,HWIN,F,DIM,WINOPT,PADOPT,WNAN)
% NDNANFILTER   N-dimensional zero-phase digital filter, ignoring NaNs.
%
%   Syntax:
%         Y = ndnanfilter(X,HWIN,F);
%         Y = ndnanfilter(X,HWIN,F,DIM);
%         Y = ndnanfilter(X,HWIN,F,DIM,WINOPT);
%         Y = ndnanfilter(X,HWIN,F,DIM,WINOPT,PADOPT);
%         Y = ndnanfilter(X,HWIN,F,DIM,WINOPT,PADOPT,WNAN);
%     [Y,W] = ndnanfilter(...);
%
%   Input:
%     X      - Data to be filtered with/without NaNs.
%     HWIN   - Window function handle (or name) or numeric multidimensional
%              window to be used (without NaNs). See WINDOW for details.
%              Default:   @rectwin  or 'rectwin' (moving average).
%     F      - A vector specifying the semi-width of the window for each
%              dimension. The final window's width will be 2*F+1.
%              Default: 3 (i.e. a 1-dimensional window of width 6).
%     DIM    - If F is a single scalar, the window will be applied through
%              this dimension; otherwise, this will be ignored.
%              Default: columns (or the first non-singleton dimension).
%     WINOPT - Cell array specifying optional arguments for the window
%              function HWIN (in addition to the width).
%              Default: {} (window's defaults).
%     PADOPT - Cell array specifying the optional arguments for the
%              PADARRAY MATLAB's function (in addition to the array X and
%              the padsize: 2*F+1). If the function is not found, data is
%              padded with zeros or the specified value: try {mean(X(:))}
%              for example.
%              Default: {'replicate'} (repeats border elements of X).
%              Default: {0} (pads with zeros if PADARRAY not found).
%     WNAN   - Integer indicating NaNs treatment and program behaviour!:
%              0: Filters data and interpolates NaNs         (default).
%              1: Filters data but do not interpolates NaNs
%              2: "Do not filters data" but interpolates NaNs!
%              See the NOTEs below
%
%   Output:
%     Y      - Filtered X data (same size as X!).
%     W      - N-dimensional window with central symmetry generated by a
%              special subfunction called NDWIND. See the description below
%              for details.
%
%   Description:
%     This function applies a N-dimensional convolution of X with W, using
%     the MATLAB's IMFILTER or CONVN function. One important aspect of the
%     function is the generation of the N-dimensional window (W) from the
%     specified function and width, which cannot be done with MATLAB's
%     functions. Besides, unlike MATLAB's FILTER, FILTER2 and IMFILTER,
%     NaNs elements are taken into account (ignored).
%
%     The N-dimensional window is generated from rotating the 1-dimensional
%     output of the HWIN function, through each of the N-dimensions, and
%     then shrinking it through each of its axes in order to fit the
%     specified semi-widths (F). This is done in the included subfunction
%     named NDWIND. In this way, the window has central symmetry and do not
%     produce a phase shift on X data.
%
%     By default, the edges are padded with the values of X at the borders
%     with the PADARRAY MATLAB's function. In this way, the edges are
%     treated smoothly. When PADARRAY is not found, the program performs
%     zero-padding.
%
%   Notes:
%     * The use of semi-widths F's is to force the generated window to be
%       even and, therefore, the change of phase is null.
%     * The window function HWIN should output an even function, otherwise,
%       it won't generate an error but the user should be aware that this
%       program will consider only the last half of it.
%     * The function window should return a monotonically decreasing
%       result, this restriction is because I try to avoid the use of FZERO
%       function, for example, to find the expanding/shrinking factors.
%     * If the user has an already generated window, it can be used in HWIN
%       instead of a function handle or name.
%     * Accepts empty value for any input. When X is empty, the program can
%       be used as a N-dimensional window generator.
%     * NaNs elements surrounded by no-NaNs elements (which will depend on
%       window width) are the ones that will be interpolated. The others
%       are leaved untouched.
%     * When WNAN=2, the programs acts like an NAN-interpolat/GAP-filling,
%       leaving untouched the no-NaNs elements but the filtering is
%       perfomed anyway. I recomend the default behaviour (WNAN=0) in order
%       to keep the filtered data in the workspace, and then use the code
%       at the end of this function to get/remove the interpolated NaNs
%       (see the example).
%     * The program looks for the IMFILTER and PADARRAY functions from the
%       Image Processing Toolbox. If not found, then CONVN is used instead
%       (slower) and pads with zeros or the given value. In this latter
%       case, if border elements are NaNs, the window won't work properly.
%
%   Example:
%     FWIN = 'hamming';
%     F = [13 8];
%     N = 100;
%     Pnoise = 0.30;
%     PNaNs  = 0.20;
%     X = peaks(N);                                     % original
%     Y = X + ((rand(size(X))-0.5)*2)*max(X(:))*Pnoise; % add noise
%     Y(round(1 + (N^2-1).*rand(N^2*PNaNs,1))) = NaN;   % add NaNs
%     [Z0,W] = ndnanfilter(Y,FWIN,F);                   % filters
%     Z1 = Z0; Z2 = Y; inan = isnan(Y);
%     Z1(inan) = NaN;
%     Z2(inan) = Z0(inan);
%     subplot(231), imagesc(X), clim = caxis; axis equal tight
%                   title('Original data')
%     subplot(232), imagesc(Y),  caxis(clim), axis equal tight
%                   title('Data + NOISE + NaNs')
%     subplot(234), imagesc(Z0), caxis(clim), axis equal tight
%                   title('FILTERS + NaNs interpolation')
%     subplot(235), imagesc(Z1), caxis(clim), axis equal tight
%                   title('FILTERS ignoring NaNs')
%     subplot(236), imagesc(Z2), caxis(clim), axis equal tight
%                   title('GAP-filling with interpolated NaNs')
%     subplot(233), imagesc(-F(1):F(1),-F(2):F(2),W), axis equal tight,
%                    title([upper(FWIN) ' 2D window']), view(2)
%
%   See also: FILTER, FILTER2 and CONVN; WINDOW from the Signal Processing
%   Toolbox; and FWIND1, FWIND2, FSPECIAL, IMFILTER and PADARRAY from the
%   Image Processing Toolbox.

%   Copyright 2008 Carlos Adrian Vargas Aguilera
%   $Revision: 1.2 $  $Date: 2008/06/30 18:00:00 $

%   Written by
%   M.S. Carlos Adrian Vargas Aguilera
%   Physical Oceanography PhD candidate
%   CICESE
%   Mexico, 2008
%   nubeobscura@hotmail.com
%
%   Download from:
%   http://www.mathworks.com/matlabcentral/fileexchange/loadAuthor.do?objec
%   tType=author&objectId=1093874

%   1.0    Release (2008/06/23 10:30:00)
%   1.1    Fixed Bug adding an extra dimension of unitary width.
%   1.2    Fixed Bug with ynan.

% Use the IMFILTER function? (faster than CONVN):
yimfilter = (exist('imfilter','file')==2);

% Use the PADARRAY function (or zero padding):
ypadarray = (exist('padarray','file')==2);

% Check inputs and sets defaults of principal arguments:
if nargin<3 || nargin>7
    error('Filtern:IncorrectNumberOfInputs',...
        'At least three inputs are needed and less than 7.')
end
if isempty(HWIN)
    HWIN = 'rectwin';
end
if isempty(F)
    F = 3;
end
N = length(F);
S = size(X);
% Secondary arguments:
if N && (nargin<4 || isempty(DIM))
    DIM = find(S~=1,1);   %  DIM = min(find(S~=1));
    if isempty(DIM), DIM = 1; end
end
if nargin<5 || isempty(WINOPT)
    WINOPT = {};
end
if nargin<6 || isempty(PADOPT)
    if ypadarray
        PADOPT = {'replicate'};
    else
        PADOPT = {0};
    end
elseif ~ypadarray && ~isnumeric(PADOPT{1})
    PADOPT = {0};
end
if nargin<7 || isempty(WNAN)
    WNAN = 0;
end

% Selects the 1-dimensional filter or set a row vector:
if N==1
    a = zeros(1,DIM);
    a(DIM) = F;
    F = a;
    clear a
end

% Checks if the window input is a function or an array:
if ~isa(HWIN,'function_handle') && ~ischar(HWIN)
    W = HWIN;
else
    W = [];
end

% If no input data but two outputs then generates the window only:
if isempty(X)
    Y = [];
    if nargout==2 && ~isempty(W)
        W = ndwind(HWIN,F,WINOPT{:});
    end
    return
end

% Generates the window:
if isempty(W)
    W = ndwind(HWIN,F,WINOPT{:});
end

% Check for NaN's:
inan = isnan(X);
ynan = any(inan(:));                       % Bug fixed 30/jun/2008
if ynan
    X(inan) = 0;
else
    factor = sum(W(:));
end

% Filtering:
if yimfilter                                % Use IMFILTER (faster)
    if ~isfloat(X)
        X = double(X);
    end
    if ~isfloat(W)
        W = double(W);
    end
    if ynan
        Y = imfilter(X,W       ,PADOPT{:},'conv');
    else
        Y = imfilter(X,W/factor,PADOPT{:},'conv');
    end
else                                        % Use CONVN
    % Sets F and S of equal sizes.
    F = reshape(F,1,N);
    Nx = numel(S);
    if N<Nx
        F(N+1:Nx) = 0;
    elseif N>Nx
        S(Nx+1:N) = 1;
    end
    F2 = 2*F;
    % Pads the borders:
    if ypadarray
        ind    = padarray(false(S),F2,true     );    % Index of the padding.
        Y      = padarray(X       ,F2,PADOPT{:});
    elseif length(PADOPT{1})==1
        ind2 = cell(N,1);
        for n = 1:N
            ind2{n} = F2(n) + (1:S(n)).';
        end
        ind          = repmat(true     ,2*F2+S);
        Y            = repmat(PADOPT{1},2*F2+S);
        ind(ind2{:}) = false;
        Y(ind2{:}) = X;
    else % No padding at all
        Y    = X;
        ind  = repmat(false,S);
        warning('Ndnanfilter:PaddingOption','Do not perfom any padding.')
    end
    % Convolutes both arrays:
    if ynan
        Y = convn(Y,W       ,'same');
    else
        Y = convn(Y,W/factor,'same');
    end
    %  Eliminates the padding:
    Y(ind) = [];
    Y      = reshape(Y,S);
end

% Estimates the averages when NaNs are present:
if ynan
    if yimfilter
        factor       = imfilter(double(~inan),W,PADOPT{:},'conv');
    else
        if ypadarray
            factor      = padarray(~inan,F2,PADOPT{:});
        elseif length(PADOPT{1})==1 % (won't work properly with NaNs at borders)
            factor          = ind;
            factor(ind2{:}) = ~inan;
        else
            factor = ~inan;
        end
        factor      = convn(factor,W,'same');
        factor(ind) = [];
        factor      = reshape(factor,S);
    end
    Y = Y./factor;
end

% What about NaNs?:
if     WNAN == 1       % Leave NaNs elements untouched!
    Y(inan) = NaN;
elseif WNAN == 2       % Leave no-NaNs elements untouched!!!
    X(inan) = Y(inan);
    Y = X;
end
end


function W = ndwind(HWIN,F,varargin)
% NDWIND Generate a N-Dimensional zero-phase window.
%
%   Syntax:
%     W = ndwind(HWIN,F);
%     W = ndwind(HWIN,F,OPT);
%
%   Input:
%     HWIN - Window function handle. See WINDOW for details. By default
%            uses: @rectwin (a rectangular window).
%     F    - A vector specifying the semiwidth of the window for each
%            dimension. The window's width will be 2*F+1. By default uses:
%            3 (i.e. a window of width 6).
%     OPT  - Cell array specifying optional arguments for the window
%            function. By default uses: {[]} (window's defaults).
%
%   Output:
%     W    - N-Dimensional window with central symmetry.
%
%   Description:
%     In the axes of each dimension, W has a 1-D window defined as
%              feval(HWIN,2*F(n)+1), n = 1,...,N.
%     That is, they are defined by the same window function but have
%     different widths. So, this program creates another widther window (at
%     least 201 points), with the same definition, and finds how much the
%     former windows should be expanded in order to fit the latter one.
%
%     Afterwards, the coordinates of every point are expanded accordingly
%     and the value of the window in those points are found by linear
%     interpolation with the bigger window.
%
%     In resume, it is like rotating this big window through every
%     dimension and then shrinking it through each of its axes to fix the
%     specified widths.
%
%   Notes:
%     * Because of the use of the semi-widths F's, all the generated
%       windows are even. Therefore the change of phase is null.
%     * The window function HWIN should output an even function, otherwise,
%       it won't generate an error but this program will consider only the
%       last half of it.
%     * The window should be monotonically decreasing.
%     * Instead of the handle window, it can be given as a string:
%       'hamming' instead of @hamming, for example.
%     * Uses the MATLAB's function FUNC2STR.
%
%   Example:
%     W = ndwind(@hamming,[3 2])
%     % Results:
%     W =
%
%              0         0    0.0800         0         0
%              0    0.1417    0.3100    0.1417         0
%              0    0.3966    0.7700    0.3966         0
%         0.0800    0.5400    1.0000    0.5400    0.0800
%              0    0.3966    0.7700    0.3966         0
%              0    0.1417    0.3100    0.1417         0
%              0         0    0.0800         0         0
%
%
%   See also: WINDOW from the Signal Processing Toolbox; and FWIND1,
%   FWIND2, and FSPECIAL from the Image Processing Toolbox.

%   Copyright 2008 Carlos Adrian Vargas Aguilera
%   $Revision: 1.1 $  $Date: 2008/06/26 19:30:00 $

%   Written by
%   M.S. Carlos Adrian Vargas Aguilera
%   Physical Oceanography PhD candidate
%   CICESE
%   Mexico, 2008
%   nubeobscura@hotmail.com
%
%   Download from:
%   http://www.mathworks.com/matlabcentral/fileexchange/loadAuthor.do?objec
%   tType=author&objectId=1093874

%   1.0    Release (2008/06/23 10:30:00)
%   1.1    Fixed Bug adding an extra dimension of unitary width.

% Check inputs:
if nargin<1 || isempty(HWIN)
    HWIN = 'rectwin';
end
if nargin<2 || isempty(F)
    F = 3;
end

% Rectangular wind?:
if isa(HWIN,'function_handle')
    HWIN = func2str(HWIN);
end
if strcmpi(HWIN,'rectwin')
    W = ones([2*F(:).'+1 1]);
    return
end

% Generate the BIG window (only the last half):
FBIG         = max([100; F(:)]);
BIGw         = feval(HWIN,2*FBIG+1,varargin{:});
BIGw(1:FBIG) = [];       % Deletes the first half.
rBIGw        = 0:FBIG;   % Window argument (distance).

% Axial windows widths:
N  = numel(F);
F  = reshape(F,1,N);
F  = [F 0];             % BUG fixed by adding an extra dimension.
N  = N+1;
F2 = 2*F+1;


% Pre-allocates the final window and the expanded axis:
W  = zeros(F2);
An = cell(N,1);
Ae = An;

% Generates the index and expanded axes:
for n = 1:N
    
    % Generate temporally the window in the n-axis:
    wn = feval(HWIN,F2(n),varargin{:});
    
    % Finds the expansion factors (Note: the window should tends to zero):
    if F(n)
        piv = wn(end);
        ind = (BIGw == piv);
        if ~any(ind)
            ind1 = (BIGw >= piv); ind1 = length(ind1(ind1));
            ind2 = (BIGw <= piv); ind2 = length(ind2(~ind2))+1;
            if ind2>FBIG+1
                r = rBIGw(ind1);
            else
                r = interp1(BIGw([ind1 ind2]), rBIGw([ind1 ind2]),piv);
            end
        else
            r = rBIGw(ind);
        end
        Ef = r/F(n);
    else
        Ef = 1;
    end
    
    % Reversed index and expanded n-axis (for the following grid):
    An{n} = (F(n):-1:0);
    Ae{n} = An{n}*Ef;
    
end

% Estimates the expanded distances outside the axes (only at the 1st
% quarter):
% Note: In a 2-Dimensional matrix, by the 1st quarter of a matrix I mean
% the first 1/4 piece of the matrix after you divided it throuh the middle
% row and column. In N-dimensions it would be the 1st 1/2^N part.
gride4      = cell(N,1);
[gride4{:}] = ndgrid(Ae{:});
R4          = sqrt(sum(reshape([gride4{:}],prod(F+1),N).^2,2));

% Generates the window and linear index in the 1st quarter:
grid4     = cell(N,1);
[grid4{:}]= ndgrid(An{:});
in        = (R4<=rBIGw(end));           % Looks for elements inside window.
W4        = zeros(F+1);                 % 1st quarter of the window.
W4(in)    = interp1(rBIGw,BIGw,R4(in)); % Interpolates the window values.
for n=1:N                               % Linear index on the 1st quarter.
    grid4{n} = flipdim(grid4{n}+1,n);
end
ind4      = sub2ind(F2,grid4{:});

% Index of permutations to fill the N-D window:
np = 2^N-1;
ip = zeros(1,np);
for n = 1:N
    ini  = 2^(n-1);
    step = ini*2;
    ip(ini:step:np) = n;
end

% Fills the N-D window by flipping W4 and the index:
ones4       = repmat(false,F2);    % Avoids using new FALSE function
ones4(ind4) = true;
W(ones4)    = W4;
for kp = ip
    W4         = flipdim(W4,kp);
    ones4      = flipdim(ones4,kp);
    W(ones4)   = W4;
end
end
