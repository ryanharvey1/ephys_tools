classdef ALLFUNCS
    %{
    Lists of included functions
    
    datastruct
    parsewilberevents
    manual_trackerjumps
    restrictMovement
    samplerate
    FixPos
    ParameterizeToLinearTrack2
    findOnsetsAndOffsets
    ndnanfilter
    ndwind
    fixNLXangle
    RemoveJumpsAndSmooth
    LoadSpikes
    ReadHeader
    FixTime
    thetamodulation
    bindata
    get_spkmatrix
    nanconv
    handle_LFP
    ThetaFilter
    RightVsLeft
    findfield
    binlaps
    FindLapsNSMAadapted
    NSMAFindGoodLaps
    peakdetz
    getPlaceFields
    kempter_lincirc
    anglereg
    infocontent
    split_trials
    split_trials_RL
    
    Ryan Harvey 2019
    %}
    methods(Static)
        function data=datastruct(path)
            ratID=strsplit(path,filesep);
            ratID=strsplit(ratID{end},'_');
            data.rat=ratID{1};
            data.sessionID=['S',strjoin(regexp(ratID{end},'\d*','Match'),'')];
            data.date_processed=date;
            data.session_path=path;
        end
        
        function out=parsewilberevents
            if exist('Events.txt','file')
                delimiter = ',';
                formatSpec = '%*q%*q%*q%f%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%q%[^\n\r]';
                fileID = fopen('Events.txt','r');
                dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
                fclose(fileID);
                eventts = dataArray{:, 1};
                eventstring = dataArray{:, 2};
                startoftask=find(ismember(eventstring,'a'));
                endoftask=find(ismember(eventstring,'b'));
                out=[eventts(startoftask(2));eventts(endoftask(2))];
            else
                out=[];
            end
        end
        
        function manual_trackerjumps(ts,x,y,StartofRec,EndofRec,path)
            % manual_trackerjumps: Allows you to manually cut out xy coordinates that
            % are outside the bounds of your maze. These can be caused by unplugs or
            % if the rat jumps out. If you do not remove these points, your ratemap
            % will be messed up.
            %
            % Input:
            %       ts
            %       x
            %       y
            %       StartofRec: ts indicating the start points of your event
            %       EndofRec: ts indicating the end points of your event
            %
            % This function depends on restrictMovement.m
            %
            % Ryan E Harvey (2018)
            
            % 1 for elminating the points outside drawn shape & 0 for inside
            restic_dir=1;
            
            savets=[];
            for i=1:length(StartofRec)
                % index out each event
                xtemp=x(ts>=StartofRec(i) & ts<=EndofRec(i));
                ytemp=y(ts>=StartofRec(i) & ts<=EndofRec(i));
                tstemp=ts(ts>=StartofRec(i) & ts<=EndofRec(i));
                % use the gui to cut out points
                [~,~,in]=ALLFUNCS.restrictMovement(xtemp,ytemp,restic_dir);
                % save the ts where the tracker error exists
                savets=[savets,tstemp(in)];
            end
            % locate the index for each tracker error
            in=ismember(ts,savets);
            % save that index to your session folder so you won't have to do this again
            % each time you run your data
            save([path,filesep,'restrictxy.mat'],'in')
        end
        
        function [x,y,in]=restrictMovement(x,y,direction)
            % restrictMovement allows you to draw a line around xy coordinates in order
            % to eliminate certain points you don't want...ie when the rat jumps out of
            % maze or tracker errors.
            %
            % Also, I have included a direction input argument so you have either
            % restrict outside or inside points. This is valuable if you have a maze
            % like a circular track where the rat could jump out away or towards the
            % center of the maze.
            %
            % Input         x,y: coordinates
            %         direction: 1 (default) to remove outside points; 0 to remove inside points
            %
            %
            % Output        x,y: retained coordinates
            %                in: logical of which coordinates were retained (so you can index ts etc.)
            %
            %
            % Ryan Harvey
            
            % check inputs
            if nargin<3
                direction=1;
            end
            
            % set up figure
            fig=figure; fig.Color=[1 1 1];plot(x,y,'k-');hold on
            title('Click around the points you want to keep')
            xlabel('X')
            ylabel('Y')
            axis equal
            disp('PRESS "ENTER" TO EXIT')
            i=1;
            
            % let the user click around the coordinates
            while true
                [X,Y]=ginput(1);
                if isempty(X)
                    break
                end
                corners(i,:)=[X,Y];
                plot(corners(:,1),corners(:,2),'r',X,Y,'*r')
                i=i+1;
            end
            
            % remove points outside or inside the shape
            in=inpolygon(x,y,corners(:,1),corners(:,2));
            if direction==1
                x=x(in);
                y=y(in);
            else
                x=x(~in);
                y=y(~in);
                in=~in;
            end
            close
        end
        
        function sampleRate=samplerate(ts)
            tempts=(ts-ts(1))/10^6;
            idx=find(tempts>=1);
            sampleRate=idx(1)-1;
        end
        
        function [x,y] = FixPos(x,y,ts,max_allowed_flips)
            % Fixes large jumps in position and interpolates missing values
            %
            % Takes all 0,0s and large jumps in position (greater than
            % jitter_threshold=15 pixels/sample) that persist for less than
            % max_allowed_flips (default = 5) and linearly interpolates the missing
            % data. Smooths conservatively afterward, as well (convolution with a gaussian, standard
            % deviation = 2 samples).
            %
            
            warning('off', 'MATLAB:interp1:NaNinY');
            
            jitter_threshold=22;
            
            bads = (x==0 | y==0);
            
            x(bads) = NaN;
            y(bads) = NaN;
            ts_=ts;
            
            if ~exist('max_allowed_flips', 'var')
                max_allowed_flips = 5; % samples
            end
            
            % [start,ends,~]=findgroups(isnan(x));
            % max_allowed_flips=max(ends-start)+1;
            
            flips = ALLFUNCS.findOnsetsAndOffsets(isnan(x));
            
            flips(:,2) = flips(:,2)+1;
            
            flips = cat(1, 1, flips(:), find([0; sqrt(diff(x).^2 + diff(y).^2)]>jitter_threshold), length(x));
            
            flips = sort(unique(flips)); % indeces of NaNs or jumps
            
            flips = [flips(1:end-1), flips(2:end)];  % epochs formation
            
            flips(:,2) = flips(:,2)-1; % adjust for diff shift
            
            flips(flips(:,2)-flips(:,1)>max_allowed_flips-1,:) = [];
            
            flips = mat2cell(flips, ones(size(flips,1),1),2); % convert to pairs corresponding to steps
            
            flips = cellfun(@(c) c(1):c(2), flips, 'unif', 0); % convert to indices
            
            x([flips{:}]) = []; % remove samples in ts and x
            y([flips{:}]) = []; % remove samples in ts and x
            ts([flips{:}]) = [];
            
            x = interp1(ts, x, ts_);
            y = interp1(ts, y, ts_);
            
            x = ALLFUNCS.ndnanfilter(x, normpdf(-6:6, 0, 2)', [], 1, {}, {}, 1); % conv with gaussian and ignore NaNs
            y = ALLFUNCS.ndnanfilter(y, normpdf(-3:3, 0, 2)', [], 1, {}, {}, 1);
            
            % fill in remaining pesky NaN values
            % XY=fillmissing([x,y],'pchip',1,'EndValues','nearest');
            % x=XY(:,1);
            % y=XY(:,2);
        end
        
        function [V,W] = ParameterizeToLinearTrack2(X,Y)
            %
            % [V,W] = ParameterizeToLinearTrack(X,Y)
            %
            % INPUTS:
            %    X,Y -- ctsd or tsd or straight arrays of data
            %
            % OUTPUTS:
            %    V,W -- ctsd or tsd or straight arrays of data
            %    V = position along track, W = orthogonal position
            %    output is of same type as input
            %
            % ALGO:
            %    Uses SVD to determine linearization.
            
            % ADR 1998
            %  version L4.0
            %  status: PROMOTED
            %
            % Updated by Ryan H 6/1/17 (lines 69-72) to stop flipping of data
            %
            %--------------------
            % Check inputs
            %--------------------
            
            if (class(X) ~= class(Y))
                error('Both or none must be tsArrays.');
            end
            if (isa(X, 'tsd') | isa(X, 'ctsd'))
                CheckTS(X, Y);
                T0 = min(StartTime(X), StartTime(Y));
                T1 = max(EndTime(X), EndTime(Y));
                X = Restrict(X,T0,T1);
                x = Data(X);
                Y = Restrict(Y,T0,T1);
                y = Data(Y);
            else
                x = X;
                y = Y;
            end
            
            %--------------------
            % Find placement of nan
            %--------------------
            
            ts=1:length(x);
            tsnan=ts(isnan(x));
            ts_not_nan=ts(~isnan(x));
            
            
            %--------------------
            % Determine eignevectors
            %--------------------
            
            Z = [x - nanmean(x)  y - nanmean(y)];
            Z=Z(~isnan(Z(:,1)),:);
            
            cv = [Z' * Z];
            [v,d] = eig(cv);
            
            A = Z * v(:,1);
            B = Z * v(:,2);
            
            %--------------------
            % Which is the long trajectory?
            %--------------------
            
            if (d(1,1) > d(2,2))
                V0 = A;
                W0 = B;
                veig=v(:,1);
                weig=v(:,2);
            else
                V0 = B;
                W0 = A;
                veig=v(:,2);
                weig=v(:,1);
            end
            
            % if (max(V0) < abs(min(V0)))
            %   % then we are reversed
            %   V0 = -V0;
            % end
            
            % zero it
            V0 = V0 - min(V0);
            
            V = V0;
            W = W0;
            
            % deal with nan values
            V=[V,ts_not_nan'];
            W=[W,ts_not_nan'];
            
            V=[V;nan(length(tsnan),1),tsnan'];
            W=[W;nan(length(tsnan),1),tsnan'];
            
            [~,I]=sort(V(:,2));
            V=V(I,1);
            W=W(I,1);
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
                    W = ALLFUNCS.ndwind(HWIN,F,WINOPT{:});
                end
                return
            end
            
            % Generates the window:
            if isempty(W)
                W = ALLFUNCS.ndwind(HWIN,F,WINOPT{:});
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
            thetaout=ALLFUNCS.RemoveJumpsAndSmooth(thetaout,max_allowed_flips,.5*pi);
            % end
            
            % fix the rest of NaN values
            % thetaout=fillmissing(thetaout,'pchip',2,'EndValues','nearest');
            
            thetaout=atan2(sin(thetaout),cos(thetaout)); % bring back to radians bound by [-pi, pi)
            
            thetaout=(thetaout*180)./pi;
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
            
            flips = ALLFUNCS.findOnsetsAndOffsets(isnan(x));
            
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
            
            % gaussian_filter = pdf('Normal',-5*smooth_std:5*smooth_std,0,smooth_std);
            if smooth_std>0
                x = ALLFUNCS.ndnanfilter(x, normpdf(filter_def, 0, smooth_std)', [], 1, {}, {}, 1); % convolve with gaussian and ignore NaNs
            end
            % x = conv(x,gaussian_filter, 'same'); % smooth by convolution with a gaussian
        end
        
        function Spikes = LoadSpikes(tfilelist)
            
            % LoadSpikes  Creates cell array of ts objects of spike times from list of tfiles
            %
            % S = LoadSpikes(tfilelist, quiet)
            %
            % INPUTS:
            %   tfilelist = cell array of strings, each of which is a tfile to open
            %       (incompatible with version unix3.1)
            %   quiet (optional):  If quiet = 1, no messages are printed
            % OUTPUTS:
            %   S = cell array such that each cell contains a ts object
            %       (timestamps which correspond to times at which the cell fired)
            %
            % ADR 1998, version L4.0, last modified '98 by ADR
            
            % status: PROMOTED
            
            %-------------------
            % Check input type
            %-------------------
            
            if ~isa(tfilelist, 'cell')
                %error('LoadSpikes: tfilelist should be a cell array.');
                tfilelist = {tfilelist};
            end
            
            nFiles = length(tfilelist);
            
            %--------------------
            % Read files
            %--------------------
            
            
            % for each tfile
            % first read the header, the read a tfile
            % note: uses the bigendian modifier to ensure correct read format.
            
            S = cell(nFiles, 1);
            for iF = 1:nFiles
                %   if ~quiet, DisplayProgress(iF, nFiles, 'Title', 'LoadSpikes'); end % REMOVED BY RYAN H 10/27/16
                tfn = tfilelist{iF};
                if ~isempty(tfn)
                    tfp = fopen(tfn, 'rb','b');
                    if (tfp == -1)
                        warning([ 'Could not open tfile ' tfn]);
                    end
                    
                    ALLFUNCS.ReadHeader(tfp);
                    S{iF} = fread(tfp,inf,'uint32');	%read as 32 bit ints
                    S{iF} = ts(S{iF});
                    
                    fclose(tfp);
                    
                end 		% if tfn valid
            end		% for all files
            for i=1:length(S)
                t=Data(S{i});
                Spikes{i}=100*t;
            end
        end
        
        function H = ReadHeader(fp)
            
            % ReadHeader  Reads NSMA header, leaves file-read-location at end of header
            %
            % H = ReadHeader(fp)
            %
            %  INPUT:
            %      fp = file-pointer to file containing header (i.e. not filename)
            %  OUTPUT:
            %      H = cell array, where each cell contains one line from the header
            %
            % ADR 1997, version L4.1, last modified '97 by ADR
            
            % status: PROMOTED
            % v4.1 17 nov 98 now works for files with no header
            
            %---------------
            % Get keys
            beginheader = '%%BEGINHEADER';
            endheader = '%%ENDHEADER';
            iH = 1; H = {};
            curfpos = ftell(fp);
            %--------------
            % go
            % look for beginheader
            headerLine = fgetl(fp);
            if strcmp(headerLine, beginheader)
                H{1} = headerLine;
                while ~feof(fp) & ~strcmp(headerLine, endheader)
                    headerLine = fgetl(fp);
                    iH = iH+1;
                    H{iH} = headerLine;
                end
            else % no header
                fseek(fp, curfpos, 'bof');
            end
        end
        
        function data=FixTime(data)
            % Aligns the session to start at time 0 & convert ts from microseconds to
            % seconds
            %
            % spikes
            for i=1:numel(data.Spikes)
                data.Spikes{i}=data.Spikes{i}-data.frames(1,1);
                data.Spikes{i}(data.Spikes{i}<0)=[];
                % convert microseconds to seconds
                data.Spikes{i}=data.Spikes{i}./10^6;
            end
            
            % lfp
            if isfield(data,'lfp') && ~isempty(data.lfp)
                data.lfp.ts=data.lfp.ts-data.frames(1,1);
                data.lfp.signal(:,data.lfp.ts<0)=[];
                data.lfp.theta(:,data.lfp.ts<0)=[];
                data.lfp.theta_phase(:,data.lfp.ts<0)=[];
                data.lfp.theta_amp(:,data.lfp.ts<0)=[];
                data.lfp.ts(data.lfp.ts<0)=[];
                % convert microseconds to seconds
                data.lfp.ts=data.lfp.ts./10^6;
            end
            % events
            if ~isempty(data.events)
                data.events=(data.events-data.frames(1,1))./10^6;
            end
            
            % trial time stamps
            data.trial_ts=(data.trial_ts-data.frames(1,1))/10^6;
            
            % uncorrected frames from linear track
            if isfield(data,'linear_track')
                data.linear_track.nonlinearFrames(:,1)=...
                    (data.linear_track.nonlinearFrames(:,1)-data.frames(1,1))./10^6;
            end
            
            %  frames
            data.frames(:,1)=(data.frames(:,1)-data.frames(1,1))./10^6;
            
            % update structure to reflect change
            data.ts_timescale='sec';
        end
        
        function [thetaindex,peak,cor,lag] = thetamodulation(spk)
            %  Calculates thetaindex per Ulanovsky
            %
            % [ind,peak,cor,lag] = b_thetaIndex(self,cel,varargin)
            
            % Constructed from the methods from
            % Grid cells without theta oscillations in the entorhinal cortex of bats
            % Michael M. Yartsev,	 Menno P. Witter	 & Nachum Ulanovsky
            % Nature 479, 103�107 (03 November 2011) doi:10.1038/nature10583
            
            % From Jason Climer and modified by Ryan E Harvey 2017
            
            ts=spk;
            
            max_lag = 0.5;
            % t_bin=0.005;
            t_bin=0.01;
            
            
            % Acor - taken from intrinsic frequency 2
            if t_bin / mod(max_lag, t_bin) ~= 2 % set lags so it is 'even' (odd number of coefficients and zero centered')
                max_lag = t_bin*floor(max_lag/t_bin)+.5*t_bin;
            end
            
            [cor, lag] = ALLFUNCS.CrossCorr(ts, ts, 'lag', [-max_lag max_lag], 'binsize', t_bin, 'norm', 'prob');
            
            % The peak of the autocorrelogram at zero lag was equalized to the maximal
            % value not including the zero lag peak
            
            
            % PERSONAL CORRESTPONANCE
            % Regarding the ylimits of the autocorrelograms - we followed the
            % convention used in some of Moshe Abeles's papers, where Abeles used
            % autocorrelograms of spike trains (Abeles, BTW, introduced the notion of
            % autocorrelations to Neuroscience in the 1970's together with George
            % Gerstein - and we all know Abeles well here in Israel... :)  ) -- and,
            % following Abeles, we normalized the min and max of the autocorrelograms
            % to be between -1 and +1
            if length(lag)~=length(cor)
                lag=linspace(-.5,.5,length(cor));
            end
            
            cor = cor/max(cor(lag~=0))*2-1;
            if any(cor>1)
                cor=rescale(cor,-1,1);
            end
            % The power spectrum of the temporal autocorrelograms was
            % assessed by computing the fast Fourier transform (FFT) of the
            % autocorrelogram, and calculating the square of the FFT magnitude
            S = abs(fft(cor)).^2;
            df = 1/range(lag); fNQ = 1/mode(diff(lag))/2;
            f = 0:df:fNQ;
            S = S(1:length(f));
            % The power spectrum was smoothed with a 2-Hz rectangular window
            S = smooth(S,2/df);
            
            % and the peak value in the 4-12 Hz band was identified
            peak = f(f>=4&f<=12);
            [~,i] = max(S(f>=4&f<=12));
            peak = peak(i);
            
            % A neuron was defined as theta)modulated if the mean power within 1)Hz of
            % each side of the peak in the 5�11 Hz frequency range was at least 5 times
            % greater than the mean spectral power between 0 Hz and 50 Hz
            thetaindex = mean(S(abs(f-peak)<1))/mean(S(f<50));
            
        end
        
        function [cor, lags] = CrossCorr(ts1, varargin)
            % [acor, lag] = CMBHOME.Spike.CrossCorr(ts1)
            % [xcor, lag] = CMBHOME.Spike.CrossCorr(ts1, ts2, varargin)
            %
            % Calculates cross correlation (or auto correlation) for vectors ts1 and
            % ts2 of spike times. ts1 and ts2 are not binned, or binary, spike trains but
            % rather the time that spikes occured.
            %
            % If ts1 is the only argument, the autocorrelation is produced (also if
            % ts1==ts2). In this case, all zeros in the set of latencies will be removed
            % before calculating the xcorr histogram.
            %
            % ARGUMENTS
            %   ts1         vector of spike times (seconds)
            %   ts2         (optional) vector of spike times
            %
            % PROPERTIES
            %   'binsize'       (default=.01 seconds) binsize for cross correlation
            %   'norm'          (default='prob') determines normalization of cross
            %                   correlogram. Can be 'prob', 'count', 'unbiased'. To use
            %                   the unbiased option, 'epoch' property must be passed.
            %   'lag'           (default [-.4 4] seconds) included to speed up algorithm.
            %                   this defines the upper and lower limit in the peristumulus
            %                   time histogram
            %   'suppress_plot' (default 1) 1 with plot, 0 without plot
            %   'epoch'         2 element vector indicating start and stop times of recording window
            %
            % RETURNS
            %   cor         a col. vector of cross correlation values corresponding to 'lag'
            %   lag         a col. vector of center-aligned lag values (seconds)
            %
            % alex wiltschko and andrew bogaard
            % Updated 3rd August 2012, Ehren Newman, increased speed for large
            % sessions.
            
            p = inputParser;
            
            p.addRequired('ts1');
            p.addOptional('ts2', [], @(c) isnumeric(c));
            p.addParamValue('norm', 'prob', @(c) ischar(c));
            p.addParamValue('binsize', .01, @(c) numel(c)==1 && (c>0));
            p.addParamValue('lag', [-.4 .4], @(c) numel(c)==2 && diff(c)>0);
            p.addParamValue('suppress_plot', 1, @(c) numel(c)==1);
            p.addParamValue('epoch', [], @(c) numel(c)==2);
            
            p.parse(ts1, varargin{:});
            
            ts1 = p.Results.ts1(:); % make col vectors
            ts2 = p.Results.ts2(:);
            norm = p.Results.norm;
            binsize = p.Results.binsize;
            lag = p.Results.lag;
            suppress_plot = p.Results.suppress_plot;
            epoch = p.Results.epoch;
            
            ac = 0;
            
            if isempty(ts2), ts2 = ts1; end
            
            if isequal(ts1,ts2), ac = 1; end % is autocorr
            
            db = nan(length(ts1), 3);
            
            s1 = 1;
            
            spkind = 1;
            
            while spkind <= length(ts1)
                
                s = s1;
                
                while ts2(s) < ts1(spkind)+lag(1) && s < length(ts2)
                    
                    s = s+1;
                    
                end
                
                s1 = s;
                
                %    f = s;
                
                f = find(ts2 <= ts1(spkind)+lag(2)==0,1);
                
                if isempty(f), f = length(ts2)+1; end
                
                %     while ts2(f) <= ts1(spkind)+lag(2)
                %
                %         f = f+1;
                %
                %         if f>length(ts2), break; end
                %
                %     end
                
                if ts2(s)<=ts1(spkind)+lag(2)
                    db(spkind, :) = [s f-1 ts1(spkind)];
                end
                
                spkind = spkind+1;
                
            end
            
            dspk = diff(db(:,1:2), 1, 2);
            
            N = 0;
            for i = 0:max(dspk)
                
                N = N + sum(dspk>=i);
                
            end
            
            
            lags = lag(1)+binsize/2:binsize:lag(2)-binsize/2;
            
            saveMem = N*4 > 1000000000;
            if saveMem
                lags = sort([lags 0-1e-10 0 0+1e-10]);
                cor = single(nan(max(dspk)+1,length(lags)));
            else
                psth = single(ones(N,1));
            end
            
            N = 1;
            for i = 0:max(dspk)
                
                where = dspk>=i;
                
                tf = db(where,1)+i;
                
                if saveMem cor(i+1,:) = single(hist(ts2(tf)-db(where,3),lags));
                else psth(N:N+sum(where)-1) = single(ts2(tf)-db(where,3)); end
                
                N = N+sum(where);
                
            end
            
            if saveMem
                cor = sum(cor);
                zeroLag = find(lags==0);
                cor = [cor(1:zeroLag-3),sum(cor(zeroLag-2:zeroLag-1),2), sum(cor(zeroLag+2:zeroLag+1),2), cor(zeroLag+3:end)];
            else
                if ac, psth(psth==0) = []; end % remove zeros in autocorrelation
                psth(psth<lags(1)-binsize/2) = [];
                cor = hist(psth, lags);
            end
            
            if strcmp(norm, 'unbiased')
                
                if isempty(epoch)
                    
                    disp('''epoch'' must be defined for unbiased normalization, no normalization performed')
                    
                else
                    
                    L = floor(diff(epoch)/binsize);
                    
                    normc = abs(lags)/binsize;
                    
                    cor = cor./(L-normc);
                    
                end
                
            elseif strcmp(norm, 'prob')
                
                cor = cor / min(length(ts1), length(ts2));
                
            end % otherwise, it will be count
            
            if ~suppress_plot
                bar(lags*1000, cor, 1, 'k'), xlim(1000*[lag(1) lag(2)]), hold on
                line([0 0], ylim, 'linestyle', ':', 'color', [.7 .7 .7]);
                set(gca, 'fontsize', 8, 'box', 'off');
            end
            
            cor = cor(:);
            
            lags = lags(:);
        end
        
        function [SmoothRateMap,nBinsx,nBinsy,occ,Coherence]=bindata(occMatrix,sampleRate,spkMatrix,track_length)
            nBinsx = round(track_length/3); nBinsy = 1;
            if isempty(occMatrix)
                SmoothRateMap=zeros(nBinsy,nBinsx);
                occ=zeros(nBinsy,nBinsx);
                Coherence=NaN;
                return
            end
            MinY = min(occMatrix(:,3));
            MaxY = max(occMatrix(:,3));
            MinX = min(occMatrix(:,2));
            MaxX = max(occMatrix(:,2));
            edges{1} = linspace(MinY, MaxY, nBinsy+1);
            edges{2} = linspace(MinX, MaxX, nBinsx+1);
            
            occMatrix = [occMatrix(:,3),occMatrix(:,2)];
            Omatrix = hist3(occMatrix,'Edges',edges);
            Omatrix(2,:) = [];
            Omatrix(:,end) = [];
            occ = Omatrix/sampleRate;
            
            % bin spike data
            Smatrix = hist3([spkMatrix(:,3), spkMatrix(:,2)],'Edges',edges);
            Smatrix(2,:) = [];
            Smatrix(:,end) = [];
            
            % divide binned spikes by occupancy to get rate maps and store data
            FilledRateMatrix = Smatrix./occ;
            FilledRateMatrix(isnan(FilledRateMatrix)) = 0;
            FilledRateMatrix(isinf(FilledRateMatrix)) = 0;
            
            % SMOOTH
            filtWidth = [1,5]; filtSigma = 1;
            imageFilter=fspecial('gaussian',filtWidth,filtSigma);
            SmoothRateMap = ALLFUNCS.nanconv(FilledRateMatrix,imageFilter, 'nanout','1d');
            
            % COHERENCE
            Coherence=corr2(FilledRateMatrix,SmoothRateMap);
        end
        
        function [spkmatrix]=get_spkmatrix(frames,ts)
            % get_spkmatrix: returns interped position values based on
            % spike times
            %
            % Ryan Harvey
            spkmatrix=[ts,interp1(frames(:,1),frames(:,2),ts),interp1(frames(:,1),frames(:,3),ts),interp1(frames(:,1),frames(:,4),ts)];
        end
        
        function c = nanconv(a, k, varargin)
            % NANCONV Convolution in 1D or 2D ignoring NaNs.
            %   C = NANCONV(A, K) convolves A and K, correcting for any NaN values
            %   in the input vector A. The result is the same size as A (as though you
            %   called 'conv' or 'conv2' with the 'same' shape).
            %
            %   C = NANCONV(A, K, 'param1', 'param2', ...) specifies one or more of the following:
            %     'edge'     - Apply edge correction to the output.
            %     'noedge'   - Do not apply edge correction to the output (default).
            %     'nanout'   - The result C should have NaNs in the same places as A.
            %     'nonanout' - The result C should have ignored NaNs removed (default).
            %                  Even with this option, C will have NaN values where the
            %                  number of consecutive NaNs is too large to ignore.
            %     '2d'       - Treat the input vectors as 2D matrices (default).
            %     '1d'       - Treat the input vectors as 1D vectors.
            %                  This option only matters if 'a' or 'k' is a row vector,
            %                  and the other is a column vector. Otherwise, this
            %                  option has no effect.
            %
            %   NANCONV works by running 'conv2' either two or three times. The first
            %   time is run on the original input signals A and K, except all the
            %   NaN values in A are replaced with zeros. The 'same' input argument is
            %   used so the output is the same size as A. The second convolution is
            %   done between a matrix the same size as A, except with zeros wherever
            %   there is a NaN value in A, and ones everywhere else. The output from
            %   the first convolution is normalized by the output from the second
            %   convolution. This corrects for missing (NaN) values in A, but it has
            %   the side effect of correcting for edge effects due to the assumption of
            %   zero padding during convolution. When the optional 'noedge' parameter
            %   is included, the convolution is run a third time, this time on a matrix
            %   of all ones the same size as A. The output from this third convolution
            %   is used to restore the edge effects. The 'noedge' parameter is enabled
            %   by default so that the output from 'nanconv' is identical to the output
            %   from 'conv2' when the input argument A has no NaN values.
            %
            % See also conv, conv2
            %
            % AUTHOR: Benjamin Kraus (bkraus@bu.edu, ben@benkraus.com)
            % Copyright (c) 2013, Benjamin Kraus
            % $Id: nanconv.m 4861 2013-05-27 03:16:22Z bkraus $
            
            % Process input arguments
            for arg = 1:nargin-2
                switch lower(varargin{arg})
                    case 'edge'; edge = true; % Apply edge correction
                    case 'noedge'; edge = false; % Do not apply edge correction
                    case {'same','full','valid'}; shape = varargin{arg}; % Specify shape
                    case 'nanout'; nanout = true; % Include original NaNs in the output.
                    case 'nonanout'; nanout = false; % Do not include NaNs in the output.
                    case {'2d','is2d'}; is1D = false; % Treat the input as 2D
                    case {'1d','is1d'}; is1D = true; % Treat the input as 1D
                end
            end
            
            % Apply default options when necessary.
            if(exist('edge','var')~=1); edge = false; end
            if(exist('nanout','var')~=1); nanout = false; end
            if(exist('is1D','var')~=1); is1D = false; end
            if(exist('shape','var')~=1); shape = 'same';
            elseif(~strcmp(shape,'same'))
                error([mfilename ':NotImplemented'],'Shape ''%s'' not implemented',shape);
            end
            
            % Get the size of 'a' for use later.
            sza = size(a);
            
            % If 1D, then convert them both to columns.
            % This modification only matters if 'a' or 'k' is a row vector, and the
            % other is a column vector. Otherwise, this argument has no effect.
            if(is1D);
                if(~isvector(a) || ~isvector(k))
                    error('MATLAB:conv:AorBNotVector','A and B must be vectors.');
                end
                a = a(:); k = k(:);
            end
            
            % Flat function for comparison.
            o = ones(size(a));
            
            % Flat function with NaNs for comparison.
            on = ones(size(a));
            
            % Find all the NaNs in the input.
            n = isnan(a);
            
            % Replace NaNs with zero, both in 'a' and 'on'.
            a(n) = 0;
            on(n) = 0;
            
            % Check that the filter does not have NaNs.
            if(any(isnan(k)));
                error([mfilename ':NaNinFilter'],'Filter (k) contains NaN values.');
            end
            
            % Calculate what a 'flat' function looks like after convolution.
            if(any(n(:)) || edge)
                flat = conv2(on,k,shape);
            else flat = o;
            end
            
            % The line above will automatically include a correction for edge effects,
            % so remove that correction if the user does not want it.
            if(any(n(:)) && ~edge); flat = flat./conv2(o,k,shape); end
            
            % Do the actual convolution
            c = conv2(a,k,shape)./flat;
            
            % If requested, replace output values with NaNs corresponding to input.
            if(nanout); c(n) = NaN; end
            
            % If 1D, convert back to the original shape.
            if(is1D && sza(1) == 1); c = c.'; end
        end
        
        function [data]=handle_LFP(data)
            % handle_LFP: loads, downsamples, and filters lfp
            %
            % ryan harvey
            
            lfpfile=dir('*.ncs');
            lfpfile={lfpfile.name};
            
            if isempty(lfpfile)
                data.lfp=[];
                return
            end
            % set up max TT number. This will make associating spike time stamps
            %   to the matching CSCs easier
            TTnum=max(str2double(extractBetween(lfpfile,'CSC','.ncs')));
            if length(TTnum)==1
                TTnum=1;
            end
            for ii=1:length(lfpfile)
                % extract CSC data with mex
                
                [Timestamps,Samples]= Nlx2MatCSC(lfpfile{ii}, [1 0 0 0 1], 0, 1);
                
                
                % reshape Samples variable into a vector
                EEGreshaped = reshape(Samples,1,length(Samples(:)));
                
                if ii==1
                    % increase the number of timestamps by a factor of 512 so that there is a timestamps for every data point
                    factor = 512;
                    EEGtsUpSampled = zeros(factor,length(Timestamps));
                    EEGtsUpSampled(1,:)=Timestamps;
                    TSdiff=diff(Timestamps);
                    TSdiff(end+1)=TSdiff(end);
                    
                    TSdiffIncrements=TSdiff/factor;
                    for i=1:length(Timestamps)
                        for j=2:factor
                            EEGtsUpSampled(j,i)=Timestamps(i)+TSdiffIncrements(i)*(j-1);
                        end
                    end
                    EEGtsUpSampled = (EEGtsUpSampled(:))';
                    
                    % downsample to 1000 Hz (is a little nicer to work with than 1523 Hz (Taube lab) or 2034 Hz (McN lab))
                    NewSFreq = 1000;
                    
                    % calculate number of samples in down sampled data
                    DurationData = (EEGtsUpSampled(end)-EEGtsUpSampled(1))/1000000; % in seconds
                    NumNSamples = floor(DurationData*NewSFreq);
                    
                    % create new timestamps
                    index = 1:NumNSamples;
                    EEG_DownSampledTimestamps(index) = EEGtsUpSampled(1)+(1000000*1/NewSFreq)*(index-1);
                    
                    % preallocate
                    EEG_DownSampledData=zeros(TTnum,length(EEG_DownSampledTimestamps));
                    EEGthetaData=zeros(TTnum,length(EEG_DownSampledTimestamps));
                    theta_phase=zeros(TTnum,length(EEG_DownSampledTimestamps));
                    theta_amp=zeros(TTnum,length(EEG_DownSampledTimestamps));
                    sec=EEG_DownSampledTimestamps/10^6;
                    sec=sec-(sec(1));
                end
                
                % check to see if ts and samples has same length (it's rare they won't)
                if length(EEGtsUpSampled)~=length(EEGreshaped)
                    [~,I]=max([length(EEGtsUpSampled);length(EEGreshaped)]);
                    if I==1
                        EEGreshaped=[EEGreshaped,zeros(1,length(EEGtsUpSampled)-length(EEGreshaped))];
                    elseif I==2
                        EEGreshaped(length(EEGtsUpSampled)+1:end)=[];
                    end
                end
                [~, filename] = fileparts(lfpfile{ii});
                trodeID = str2double(extractAfter(filename,'CSC'));
                % interpolate data points
                EEG_DownSampledData(TTnum,:)=interp1(EEGtsUpSampled, EEGreshaped, EEG_DownSampledTimestamps);
                EEGthetaData(TTnum,:) = ALLFUNCS.ThetaFilter(EEG_DownSampledData(TTnum,:),NewSFreq);
                
                % RUN FMA PHASE CODE
                [phase,amplitude,~]=Phase([sec',EEGthetaData(TTnum,:)']);
                theta_phase(TTnum,:)=phase(:,2)';
                theta_amp(TTnum,:)=amplitude(:,2)';
                %     ts_sec=phase(:,1)';
                
                clearvars -except sampleout lfpfile i EEG_DownSampledTimestamps EEGtsUpSampled...
                    EEG_DownSampledData NewSFreq EEGthetaData sec theta_phase theta_amp ts_sec data
            end
            % save([extractBefore(eegfile{1,:},'CSC'),'LFP.mat'],'EEG_DownSampledTimestamps',...
            %     'EEG_DownSampledData','EEGthetaData','theta_phase','theta_amp','-v7.3')
            
            data.lfp.ts=EEG_DownSampledTimestamps;
            data.lfp.signal=EEG_DownSampledData;
            data.lfp.theta=EEGthetaData;
            data.lfp.theta_phase=theta_phase;
            data.lfp.theta_amp=theta_amp;
        end
        
        function signal_filtered = ThetaFilter(signal, Fs)
            % Takes 'signal' and bandpasses it to theta frequencies (6 to 10 Hz)
            %
            % Arguments
            %
            % signal - arbitrary signal to be filtered
            % Fs - sample frequency
            %
            % signal_filtered = ThetaFilter(signal, Fs)
            
            Wn_theta = [4/(Fs/2) 12/(Fs/2)]; % normalized by the nyquist frequency
            
            [btheta,atheta] = butter(3,Wn_theta);
            
            signal_filtered = filtfilt(btheta,atheta,signal);
        end
        
        function [right,left,DirectionalityIndex,Displacement,nlaps,rateoverlap,...
                fieldoverlap,spatialcorrelation,startgoodrun,stopgoodrun,laps]=...
                RightVsLeft(data_video_nospk,data_video_spk,track_length,sampleRate,data)
            
            if isfield(data.linear_track,'lapinfo')
                startgoodrun=data.linear_track.lapinfo.startgoodrun;
                stopgoodrun=data.linear_track.lapinfo.stopgoodrun;
                laps=data.linear_track.lapinfo.laps;
            else
                % locate laps
                %     laps=FindLapsNSMAadapted(data_video_nospk(:,1),rescale(data_video_nospk(:,2),0,track_length),35);
                % locate good laps
                %     [startgoodrun,stopgoodrun,laps]=NSMAFindGoodLaps(data_video_nospk(:,1),rescale(data_video_nospk(:,2),0,360),laps);
                %
                % if <15 good laps exist, dynamically loosening the edge (.1 to .3) completeness (.1 .3) threshold criteria
                % This normally should not be an issue if your animal is making nice ballistic laps
                % and your tracker is accurate.
                % note that ideally your rat should make >10 laps each session
                %     if size([laps.start_ts],2)/2<15
                if track_length<100
                    posbins=40;
                else
                    posbins=50;
                end
                param=[[nchoosek((1:3),1),nchoosek((1:3),1)];nchoosek((1:3),2);nchoosek((3:-1:1),2)]*0.1;
                for i=1:length(param)
                    laps=ALLFUNCS.FindLapsNSMAadapted(data_video_nospk(:,1),rescale(data_video_nospk(:,2),0,track_length),35);
                    [~,~,laps]=ALLFUNCS.NSMAFindGoodLaps(data_video_nospk(:,1),rescale(data_video_nospk(:,2),0,track_length),laps,...
                        param(i,1),param(i,2),0,posbins);
                    nlaps(i,1)=size([laps.start_ts],2)/2;
                end
                [~,I]=max(nlaps);
                laps=ALLFUNCS.FindLapsNSMAadapted(data_video_nospk(:,1),rescale(data_video_nospk(:,2),0,track_length),35);
                [startgoodrun,stopgoodrun,laps]=ALLFUNCS.NSMAFindGoodLaps(data_video_nospk(:,1),rescale(data_video_nospk(:,2),0,track_length),laps,...
                    param(I,1),param(I,2),0,posbins);
                %     end
            end
            
            if length([laps.direction])==1
                right.occ=nan(1,track_length/3);
                right.nBinsx=track_length/3;
                right.nBinsy=1;
                right.Coherence=NaN;
                right.SmoothRateMap=nan(1,track_length/3);
                right.dataspks=nan(1,6);
                right.lap_perm_stability=NaN;
                right.stabilityoverlaps=NaN;
                right.meanstability=NaN;
                left.occ=nan(1,track_length/3);
                left.nBinsx=track_length/3;
                left.nBinsy=1;
                left.Coherence=NaN;
                left.SmoothRateMap=nan(1,track_length/3);
                left.dataspks=nan(1,6);
                left.lap_perm_stability=NaN;
                left.stabilityoverlaps=NaN;
                left.meanstability=NaN;
                DirectionalityIndex=NaN;
                Displacement=NaN;
                nlaps=1;
                rateoverlap=NaN;
                fieldoverlap=NaN;
                spatialcorrelation=NaN;
                return
            end
            
            % number of laps
            nlaps=sum([laps.direction]==laps(2).direction);
            
            % find times where rat turns around without finishing a lap to remove later
            % [startgoodrun,stopgoodrun]=FindTurnAroundsNSMAadapted(data_video_nospk,laps,2);
            % find laps in spike matrix
            
            lefts=[];rights=[];i=1;
            while true
                if laps(i).direction==1
                    lefts=[lefts;data_video_spk(find(data_video_spk(:,1)==laps(i).start_ts):...
                        find(data_video_spk(:,1)==laps(i+1).start_ts),:)];
                elseif laps(i).direction==-1
                    rights=[rights;data_video_spk(find(data_video_spk(:,1)==laps(i).start_ts):...
                        find(data_video_spk(:,1)==laps(i+1).start_ts),:)];
                end
                i=i+1;
                if i>=length(laps)
                    break
                end
            end
            % remove turn around times noted above
            tempdata=[rights;lefts];
            [~,sortI]=sort(tempdata(:,1));
            sorted=tempdata(sortI,:);
            stopgoodrun(end)=sorted(end,1);
            [~,Is]=intersect(sorted(:,1),startgoodrun);
            [~,Ie]=intersect(sorted(:,1),stopgoodrun);
            c=size(sorted,2)+1;
            for i=1:length(Is)
                sorted(Is(i):Ie(i),c)=ones(length(Is(i):Ie(i)),1);
            end
            newInd(sortI) = 1:length(tempdata);
            unsorted=sorted(newInd,:);
            rights=unsorted(1:length(rights),:);
            lefts=unsorted(length(rights)+1:end,:);
            rights(rights(:,end)==0,:)=[];
            lefts(lefts(:,end)==0,:)=[];
            rights(:,end)=[];
            lefts(:,end)=[];
            
            % put laps in right and left structs
            right.datavid=rights(rights(:,end)==0,:);
            right.dataspks=rights;
            right.spks_VEL=rights(rights(:,end)==1,:);
            right.num_spikes=size(right.spks_VEL,1);
            
            left.datavid=lefts(lefts(:,end)==0,:);
            left.dataspks=lefts;
            left.spks_VEL=lefts(lefts(:,end)==1,:);
            left.num_spikes=size(left.spks_VEL,1);
            
            % BIN DATA
            [right.SmoothRateMap,right.nBinsx,right.nBinsy,right.occ,right.Coherence]=...
                ALLFUNCS.bindata(right.datavid,sampleRate,right.spks_VEL,track_length);
            
            % SAME COMMENTS AS ABOVE, BUT OPPOSITE DIRECTION
            [left.SmoothRateMap,left.nBinsx,left.nBinsy,left.occ,left.Coherence]=...
                ALLFUNCS.bindata(left.datavid,sampleRate,left.spks_VEL,track_length);
            
            % CALCULATE DIRECTIONALITY INDEX (Ravassard, 2013)
            % DirectionalityIndex=(right.SmoothRateMap-left.SmoothRateMap)/(right.SmoothRateMap+left.SmoothRateMap);
            DirectionalityIndex=abs(sum(right.SmoothRateMap-left.SmoothRateMap)/sum(right.SmoothRateMap+left.SmoothRateMap));
            
            % LINEAR DISPLACEMENT (Dispalcement unit is in cm)
            % for ibins=0:length(right.SmoothRateMap)-1
            %     SmoothRateMap_Left=circshift(left.SmoothRateMap,[0 ibins]);
            %     Correlations(ibins+1,:)=corr2(right.SmoothRateMap,SmoothRateMap_Left);
            % end
            % [~,Displacement]=max(Correlations); Displacement=Displacement-1;
            
            % SPATIAL CORRELATION
            spatialcorrelation=corr2(right.SmoothRateMap,left.SmoothRateMap);
            
            % RATE OVERLAP
            normratemaps=rescale([right.SmoothRateMap,left.SmoothRateMap],0,1);
            norm1=normratemaps(1:length(right.SmoothRateMap));
            norm2=normratemaps(length(right.SmoothRateMap)+1:end);
            rateoverlap=abs(max(norm1)-max(norm2));
            
            % FIELD OVERLAP
            [~,peak1]=max(right.SmoothRateMap);
            [~,peak2]=max(left.SmoothRateMap);
            if peak1==peak2
                fieldoverlap=1;
            else
                [field1]=ALLFUNCS.findfield(right.SmoothRateMap);
                [field2]=ALLFUNCS.findfield(left.SmoothRateMap);
                fieldoverlap=length(intersect(find(field1),find(field2)))/sum([field1+field2]>0);
            end
            
            % DISPLACEMENT
            [~,i1]=max(right.SmoothRateMap);[~,i2]=max(left.SmoothRateMap);
            Displacement=abs(i1-i2)*(track_length/length(left.SmoothRateMap));
            
            % RIGHT
            % LAP STABILITY MEASURES
            
            i=1;
            mapsL=[];
            mapsR=[];
            corsL=[];
            corsR=[];
            lapsL=[];
            lapsR=[];
            while true
                if laps(i).direction==1
                    lapsL{i,1}=left.dataspks(find(left.dataspks(:,1)==laps(i).start_ts):find(left.dataspks(:,1)==laps(i+1).start_ts),:);
                    map=ALLFUNCS.binlaps(left.dataspks,lapsL{i},sampleRate,lapsL{i}(lapsL{i}(:,end)==1,:),track_length);
                    corsL{i,1}=corr2(left.SmoothRateMap,map);
                    mapsL{i,1}=map;
                elseif laps(i).direction==-1
                    lapsR{i,1}=right.dataspks(find(right.dataspks(:,1)==laps(i).start_ts):find(right.dataspks(:,1)==laps(i+1).start_ts),:);
                    map=ALLFUNCS.binlaps(right.dataspks,lapsR{i},sampleRate,lapsR{i}(lapsR{i}(:,end)==1,:),track_length);
                    corsR{i,1}=corr2(right.SmoothRateMap,map);
                    mapsR{i,1}=map;
                end
                i=i+1;
                if i>=length(laps)
                    if isempty(lapsL)
                        left.cors=NaN;
                        left.maps=NaN;
                        left.laps=NaN;
                    else
                        left.cors=corsL(~cellfun(@isempty,lapsL));
                        left.maps=mapsL(~cellfun(@isempty,lapsL));
                        left.laps=lapsL(~cellfun(@isempty,lapsL));
                    end
                    if isempty(lapsR)
                        right.cors=NaN;
                        right.maps=NaN;
                        right.laps=NaN;
                    else
                        right.cors=corsR(~cellfun(@isempty,lapsR));
                        right.maps=mapsR(~cellfun(@isempty,lapsR));
                        right.laps=lapsR(~cellfun(@isempty,lapsR));
                    end
                    break
                end
            end
            
            % ALL PERM CORRELATIONS OF LAPS CORRELATED TO EACHOTHER
            if isempty(lapsR)
                right.lap_perm_stability=NaN;
            else
                corrz=[];
                for m=1:size(right.maps,1)
                    for mm=1:size(right.maps,1)
                        corrz=[corrz;corr2(right.maps{m},right.maps{mm})];
                    end
                end
                right.lap_perm_stability=nanmean(corrz);
            end
            if isempty(lapsL)
                left.lap_perm_stability=NaN;
            else
                corrz=[];
                for m=1:size(left.maps,1)
                    for mm=1:size(left.maps,1)
                        corrz=[corrz;corr2(left.maps{m},left.maps{mm})];
                    end
                end
                left.lap_perm_stability=nanmean(corrz);
            end
            
            % STABILITY OVER LAPS: SEE IF LAPS ARE SIMILIARITY SIMILAR TO THE FINAL RESULTING RATE MAP OR NOT
            % + CORR: BECOME STABLE
            % - CORR: LOSE STABILITY
            % 0 CORR: ARE STABLE THROUGHOUT LAPS
            if isempty(lapsR)
                right.stabilityoverlaps=NaN;
            else
                try
                    right.stabilityoverlaps=corr2([1:length(right.cors)]',[right.cors{:}]');
                catch
                    right.stabilityoverlaps=NaN;
                end
            end
            if isempty(lapsL)
                left.stabilityoverlaps=NaN;
            else
                try
                    left.stabilityoverlaps=corr2([1:length(left.cors)]',[left.cors{:}]');
                catch
                    left.stabilityoverlaps=NaN;
                end
            end
            % MEAN STABILITY
            if isempty(lapsR)
                right.meanstability=NaN;
            else
                right.meanstability=nanmean([right.cors{:}]);
            end
            if isempty(lapsL)
                left.meanstability=NaN;
            else
                left.meanstability=nanmean([left.cors{:}]);
            end
        end
        
        function [field]=findfield(Ratemap)
            [M,I]=max(Ratemap); Thres=M*0.20;
            field=zeros(1,length(Ratemap)); field(I)=1;
            
            % FORWARD
            for i=1:length(Ratemap)-I
                if Ratemap(I)~=length(Ratemap)
                    if Ratemap(I+i)>Thres
                        field(I+i)=1;
                    elseif Ratemap(I+i)<Thres
                        break
                    end
                end
            end
            
            % BACKWARD
            for i=1:I-1
                if Ratemap(I)~=1
                    if Ratemap(I-i)>Thres
                        field(I-i)=1;
                    elseif Ratemap(I-i)<Thres
                        break
                    end
                end
            end
        end
        
        function [SmoothRateMap,occ]=binlaps(occ_overall,occMatrix,sampleRate,spks_VEL,track_length)
            
            nBinsx = round(track_length/3); nBinsy = 1;
            
            if isempty(occMatrix)
                SmoothRateMap=zeros(nBinsy,nBinsx);
                return
            end
            
            MinY = min(occ_overall(:,3));
            MaxY = max(occ_overall(:,3));
            MinX = min(occ_overall(:,2));
            MaxX = max(occ_overall(:,2));
            edges{1} = linspace(MinY, MaxY, nBinsy+1);
            edges{2} = linspace(MinX, MaxX, nBinsx+1);
            
            occMatrix = [occMatrix(:,3),occMatrix(:,2)];
            Omatrix = hist3(occMatrix,'Edges',edges);
            Omatrix(2,:) = [];
            Omatrix(:,end) = [];
            occ = Omatrix/sampleRate;
            
            % bin spike data
            Smatrix = hist3([spks_VEL(:,3), spks_VEL(:,2)],'Edges',edges);
            Smatrix(2,:) = [];
            Smatrix(:,end) = [];
            
            % divide binned spikes by occupancy to get rate maps and store data
            % (can use occ+eps instead of removing 0's)
            FilledRateMatrix = Smatrix./occ;
            FilledRateMatrix(isnan(FilledRateMatrix)) = 0;
            FilledRateMatrix(isinf(FilledRateMatrix)) = 0;
            
            % SMOOTH
            filtWidth = [1,5]; filtSigma = 1;
            imageFilter=fspecial('gaussian',filtWidth,filtSigma);
            SmoothRateMap = ALLFUNCS.nanconv(FilledRateMatrix,imageFilter, 'nanout','1d');
        end
        
        function laps=FindLapsNSMAadapted(Vts,Vdata,newLapThreshold)
            %
            % [laps, Vhs] = FindLaps_HorseShoe(V,newLapThreshold);
            %
            % Find Laps in "Horseshoe-Geometry" of a circular track
            %
            % INPUT:
            % V    ...  tsd of angular position in [0,360] degree range, with horsehoe
            %           opening at 0 degrees.
            %
            % newLapThreshold ... OPTIONAL endpoint proximity threshold in percent of track length (default = 15%);
            %                     whenever rat enters the proximity zone of e.g. 15% of tracklength near a horseshoe end, a new lap
            %                     is started and the maximum (or minimum) is searched
            %                     for a Lap-Top (around 360 end)  or Lap-Bottom (around 0 end).
            %
            % OUTPUT:
            % laps  .... 1*nLaps struct array with fields
            %   laps(i).start_ts  ... start timestamp of i-th lap (in 0.1 millisec
            %                              NSMA units)
            %   laps(i).pos       ... the value of input position V (in circle geometry not horseshoe!) at lap start point (in degrees)
            %   laps(i).start_idx ... the index of the new lap start frame in input V tsd
            %   laps(i).direction ... +1/-1 for up/down laps (with increasing/decreasing position in angle degrees)
            %
            %  Vhs ....  tsd nFrames*1 vector of position in horseshoe geometry, where
            %               negative direction laps are ranging from [-360,0] and positive direction
            %               laps are ranging from [0 360]. NOTE: there are discontinuous jumps at
            %               the horseshoe turn points around 0 degrees and around
            %               360/-360 degrees. NOTE2: nFrames = length(Data(V))
            %
            %
            %   Author: PL
            %   Version: 0.9  05/12/2005
            %   edited by Ryan Harvey to work with standard linear track
            
            if nargin == 1
                newLapThreshold = 15;   % default value
            end
            
            % Vdata = data_video(:,2);
            % Vts = data_video(:,1);
            
            TL = abs(max(Vdata) - min(Vdata));   % track length in degrees
            th1= min(Vdata)+TL*newLapThreshold/100;        % lower threshold for lower horeshoe end (degrees)
            th2 = max(Vdata)-TL*newLapThreshold/100;      % upper threshold for upper horseshoe end (degrees)
            
            
            % loop over all frames
            laps(1).start_ts = Vts(1);
            laps(1).pos = Vdata(1);
            laps(1).start_idx = 1;
            laps(1).direction = 0;
            iLap = 1;
            
            newUpThCross = 1;     % flag for new lap top search
            newDownThCross = 1;     % flag for new lap top search
            for i=1:length(Vdata)
                if Vdata(i) < th1    % search for min
                    if newUpThCross       % start a new lap
                        newUpThCross = 0;
                        newDownThCross = 1;
                        iLap = iLap + 1;
                        laps(iLap).start_ts = Vts(i);
                        laps(iLap).pos = Vdata(i);
                        laps(iLap).start_idx = i;
                        laps(iLap).direction = +1;
                    end
                    if Vdata(i) < laps(iLap).pos       % record new min if any
                        laps(iLap).start_ts = Vts(i);
                        laps(iLap).pos = Vdata(i);
                        laps(iLap).start_idx = i;
                    end
                end
                
                if Vdata(i) > th2   % search for max
                    if newDownThCross       % start a new lap
                        newUpThCross = 1;
                        newDownThCross = 0;
                        iLap = iLap + 1;
                        laps(iLap).start_ts = Vts(i);
                        laps(iLap).pos = Vdata(i);
                        laps(iLap).start_idx = i;
                        laps(iLap).direction = -1;
                    end
                    if Vdata(i) > laps(iLap).pos       % record new min if any
                        laps(iLap).start_ts = Vts(i);
                        laps(iLap).pos = Vdata(i);
                        laps(iLap).start_idx = i;
                    end
                    
                end
                
            end
            
            % fix direction of first lap which was unknown above
            % make first lap direction opposite of second lap's direction (laps alternate!)
            laps(1).direction = -laps(2).direction;
            
            % make sure laps cross the halfway point
            middle=median(min(Vdata):max(Vdata));
            i=1;
            while true
                if any(laps(i).pos:laps(i+1).pos>middle) && any(laps(i).pos:laps(i+1).pos<middle)==0
                    laps(i+1)=[];
                end
                i=i+1;
                if i>=length([laps.pos])
                    if length([laps.pos])<iLap
                        laps(1).direction = -laps(2).direction;
                    end
                    break
                end
            end
        end
        
        function [startgoodlaps,stopgoodlaps,laps]=NSMAFindGoodLaps(ts,V_rest,laps,edgethresh,completeprop,plotlaps,posbins)
            % [startgoodlaps, stopgoodlaps, laps] =
            %        FindGoodLaps(V_rest,laps,edgethresh,completeprop,plotlaps,posbins)
            %
            % find and eliminate laps which have too many NaNs (because rat was off
            % track), and parts of laps where rat turns around in middle of track
            %
            % inputs: V_rest: V coordinates of rat with off track periods masked out
            %                 (as NaNs)
            %         laps: struct with lap start and end times (generated by
            %               FindLaps_Horseshoe)
            %         edgethresh: threshold for detection of a turn around point
            %                     (proportion of length of track) (default = 0.1)
            %         completeprop: the amount of lap that can be missing (NaNs) to
            %                       still be considered a lap (default = 0.2).
            %         plotlaps: flag for making plots of each lap, and pause for user
            %                   to hit key to continue (default = 1)
            %         posbins: number of bins to divide the track into to determine
            %                  position coverage percentage; at 60frames/s want at
            %                  least 2cm/bin (default = 50bins; this works for 100+ cm
            %                  track, as long as V_rest is in cm)
            % outputs: startgoodlaps, stopgoodlaps: start and stop times of good lap
            %                                       periods
            %          laps: a new laps struct, with the bad laps removed
            %
            % ZN 04/2011
            
            if nargin==3
                edgethresh=0.1;
                completeprop=0.2;
                plotlaps=0;
                posbins=50;
            elseif nargin==4
                completeprop=0.2;
                plotlaps=1;
                posbins=50;
            elseif nargin==5
                plotlaps=1;
                posbins=50;
            elseif nargin==6
                posbins=50;
            end
            
            if edgethresh>1      % in case edgethresh is input as a percentage instead of a proportion
                edgethresh=edgethresh/100;
            end
            if completeprop>1      % in case completeprop is input as a percentage instead of a proportion
                completeprop=completeprop/100;
            end
            
            bottomend=min((V_rest));
            topend=max((V_rest));
            bins=bottomend:(topend-bottomend)/posbins:topend;
            delta=(topend-bottomend)*edgethresh;     % threshold for peak/trough detection
            startgoodlaps=[];
            stopgoodlaps=[];
            if plotlaps
                a=figure;
            end
            l=1;
            while l<=length(laps) %length(laps)
                % select out just this lap
                if l==length(laps)
                    endoflap=ts(end);
                else
                    endoflap=laps(l+1).start_ts;
                end
                %     Vlap=Restrict(V_rest, laps(l).start_ts, endoflap);
                %     t=Range(Vlap);
                %     v=Data(Vlap);
                
                v= V_rest(find(ts==laps(l).start_ts):find(ts==endoflap));
                t= ts(find(ts==laps(l).start_ts):find(ts==endoflap));
                
                
                if plotlaps
                    figure(a); subplot(1,5,1:4); plot(t,v,'k-')
                    if length(t)>1
                        axis([t(1) t(end) bottomend topend])
                    end
                    title(['lap # ', num2str(l), ': ', num2str(laps(l).direction)])
                    xlabel('Time')
                    ylabel('Position')
                end
                
                % find turn around points during this lap
                lookformax=laps(l).direction==1;
                [peak,trough]=ALLFUNCS.peakdetz(v, delta, lookformax, 0);
                if plotlaps
                    hold on;
                    if ~isempty(peak)
                        plot(t(peak(:,1)), peak(:,2), 'r.')
                    end
                    if ~isempty(trough)
                        plot(t(trough(:,1)), trough(:,2), 'b.')
                    end
                end
                
                if lookformax
                    % find the direct path from bottomend to topend (or mark lap for
                    % deleting if the turn around points are not in those ranges)
                    if ~isempty(trough)
                        % find the last trough in range of bottomend (start of lap)
                        [gt,blah]=size(trough);
                        while gt>0 && trough(gt,2)>=2*delta+bottomend
                            gt=gt-1;
                        end
                        % assign the next peak after that trough as the end of the lap
                        % (or mark lap for deleting, if that peak is not at topend)
                        if gt==0
                            if peak(1,2)>topend-2*delta
                                t=t(1:peak(1,1));
                                v=v(1:peak(1,1));
                            else
                                % this marks the lap for deleting
                                t=t(1:5);
                                v=v(1:5);
                                if plotlaps
                                    title('continous run does not reach end of track')
                                end
                            end
                        else
                            [et,blah]=size(peak);
                            if gt+1>et
                                gt=0;
                                t=t(1:2);
                                v=v(1:2);
                                if plotlaps
                                    title('last start of run is past last peak')
                                    % this happens if the lap never leaves the start of lap region
                                end
                            else
                                t=t(trough(gt,1):peak(gt+1,1));
                                v=v(trough(gt,1):peak(gt+1,1));
                            end
                        end
                        if plotlaps
                            plot(t, v, 'g--')
                        end
                        
                    else % if ~isempty(trough)
                        % make sure peak exists and is in range of topend
                        if isempty(peak)
                            if plotlaps
                                title('Peak does not exist')
                            end
                            if length(t)>2
                                t=t(1:2);
                                v=v(1:2);
                            end
                        elseif peak(1,2)<topend-2*delta
                            % this marks the lap for deleting
                            if length(t)>5
                                t=t(1:5);
                                v=v(1:5);
                            end
                            if plotlaps
                                title('run does not end at end of track ')
                            end
                        end
                    end
                else % if lookformax
                    % find the direct path from topend to bottomend (or mark lap for
                    % deleting if the turn around points are not in those ranges)
                    if ~isempty(peak)
                        % find the last peak in range of topend (start of lap)
                        [gt,blah]=size(peak);
                        while gt>0 && peak(gt,2)<=topend-2*delta
                            gt=gt-1;
                        end
                        % assign the next trough after that peak as the end of the lap
                        % (or mark lap for deleting, if that trough is not at bottomend)
                        if gt==0
                            if trough(1,2)<bottomend+2*delta
                                t=t(1:trough(1,1));
                                v=v(1:trough(1,1));
                            else
                                % this marks the lap for deleting
                                t=t(1:5);
                                v=v(1:5);
                                if plotlaps
                                    title('continous run does not reach end of track')
                                end
                            end
                        else
                            [et,blah]=size(trough);
                            if gt+1>et
                                t=t(1:2);
                                v=v(1:2);
                                gt=0;
                                if plotlaps
                                    title('last start of run is past last trough')
                                    % this happens if the lap never leaves the start of lap region
                                end
                            else
                                t=t(peak(gt,1):trough(gt+1,1));
                                v=v(peak(gt,1):trough(gt+1,1));
                            end
                        end
                        if plotlaps
                            plot(t, v, 'g--')
                        end
                    else % if ~isempty(peak)
                        % make sure trough exists and is in range of bottomend
                        if isempty(trough)
                            if plotlaps
                                title('Through does not exist')
                            end
                            if length(t)>2
                                t=t(1:2);
                                v=v(1:2);
                            end
                        elseif trough(1,2)>bottomend+2*delta
                            % this marks the lap for deleting
                            if length(t)>5
                                t=t(1:5);
                                v=v(1:5);
                            end
                            if plotlaps
                                title('run does not end at end of track ')
                            end
                        end
                    end
                end % else looformax
                
                vcovered=hist(v,bins);
                if plotlaps
                    subplot(1,5,5); plot(vcovered, bins,'k-')
                    histpeak=max([vcovered,1]);
                    axis([0 histpeak bottomend topend])
                    xlabel('frames spent in position')
                end
                
                if length(v)<3
                    % eliminate the lap if it is non-existent (as is sometimes the case for lap 1)
                    laps(l)=[];
                    if plotlaps
                        title(['Non-existent lap. ', num2str(length(laps)), ' laps left.'])
                    end
                    % eliminate the lap if >completeprop of it is NaNs or if it has been marked for
                    % deleting above
                elseif length(v)<6 || sum(vcovered==0)>completeprop*posbins
                    laps(l)=[];
                    % remove the other lap from the lap pair
                    if mod(l,2)==0
                        % delete previous lap from laps
                        laps(l-1)=[];
                        if plotlaps
                            title(['Deleting this and the previous lap. ', num2str(length(laps)), ' laps left.'])
                        end
                        % change goodlaps markers to delete previous lap from laps
                        if ~isempty(stopgoodlaps)
                            if ~exist('lastlapend', 'var') || startgoodlaps(end)>lastlapend
                                startgoodlaps(end)=[];
                                stopgoodlaps(end)=[];
                            else
                                stopgoodlaps(end)=lastlapend;
                            end
                        end
                        l=l-1;
                    elseif l<=length(laps) && l>1
                        % delete next lap from laps
                        laps(l)=[];
                        if plotlaps
                            title(['Deleting this and the next lap. ', num2str(length(laps)), ' laps left.'])
                        end
                    else
                        if plotlaps
                            title(['Deleting this lap. ', num2str(length(laps)), ' laps left.'])
                        end
                    end
                else % if lap is good
                    % store last lap end just in case have to delete this lap with next lap
                    if ~isempty(stopgoodlaps)
                        lastlapend=stopgoodlaps(end);
                    end
                    % add this lap to goodlaps
                    if ~isempty(stopgoodlaps) && stopgoodlaps(end)==t(1)
                        stopgoodlaps(end)=t(end);
                    else
                        startgoodlaps=[startgoodlaps; t(1)];
                        stopgoodlaps=[stopgoodlaps; t(end)];
                    end
                    l=l+1;
                end
                if plotlaps
                    drawnow
                    %         pause
                    hold off;
                end
            end
        end
        function [maxtab, mintab]=peakdetz(v, delta, lookformax, backwards)
            %PEAKDET Detect peaks in a vector
            %        [MAXTAB, MINTAB] = PEAKDETZ(V, DELTA, lookformax, backwards) finds
            %        the local maxima and minima ("peaks") in the vector V.
            %        A point is considered a maximum peak if it has the maximal
            %        value, and was preceded (to the left) by a value lower by
            %        DELTA. MAXTAB and MINTAB consists of two columns. Column 1
            %        contains indices in V, and column 2 the found values.
            %
            % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
            % This function is released to the public domain; Any use is allowed.
            %
            % ZN edit 04/2010: added option to specify looking for troughs or peaks
            % first (lookformax variable: if 1, will look for peaks first, if 0 will
            % look for troughs; default is look for peaks); and option to go backwards
            % (so that find last instance of a peak/trough value instead of the first
            % instance: backwards variable: if 1 will go backwards, if 0 or absent,
            % will go forwards); and changed it so that last min/max value will be
            % assigned
            
            if nargin<3
                lookformax = 1;
            end
            
            maxtab = [];
            mintab = [];
            
            v = v(:); % Just in case this wasn't a proper vector
            
            if (length(delta(:)))>1
                error('Input argument DELTA must be a scalar');
            end
            
            if delta <= 0
                error('Input argument DELTA must be positive');
            end
            
            if nargin<4 || backwards==0
                inc=1;
                first=1;
                last=length(v);
            elseif backwards
                inc=-1;
                first=length(v);
                last=1;
            end
            
            
            mn = Inf; mx = -Inf;
            mnpos = NaN; mxpos = NaN;
            
            for ii=first:inc:last
                this = v(ii);
                if this > mx, mx = this; mxpos = ii; end
                if this < mn, mn = this; mnpos = ii; end
                
                if lookformax
                    if this < mx-delta || (ii==last && ~isempty(mintab) && mx-delta>mintab(end,2))
                        maxtab = [maxtab ; mxpos mx];
                        mn = this; mnpos = ii;
                        lookformax = 0;
                    end
                else
                    if this > mn+delta || (ii==last && ~isempty(maxtab) && mn+delta<maxtab(end,2))
                        mintab = [mintab ; mnpos mn];
                        mx = this; mxpos = ii;
                        lookformax = 1;
                    end
                end
            end
            
            if isempty([mintab;maxtab])
                if lookformax
                    if mx-mn>delta
                        maxtab=[mxpos mx];
                    end
                else
                    if mx-mn>delta
                        mintab=[mnpos mn];
                    end
                end
            end
        end
        
        function [fields]=getPlaceFields(varargin)
            % USAGE
            %
            %
            % INPUTS
            %
            %   ratemap             MxN matrix where M is the number of cells, N is the number
            %                       of spatial bins
            %   minPeakRate         minimum rate for peak of field [default: 2]
            %   minFieldWidth       minimum width of field [default: 2]
            %   maxFieldWidth       maximum width of field [default: 30]
            %   percentThreshold    percent change between peak rate and start/stop of field
            %
            %
            % OUTPUTS
            %
            %   fields struct with field data for each cell
            %
            %
            %
            % HELP
            % This function tries to identify place fields based on a set of heuristics
            %
            % written by david tingley, 2017
            % Adapted by Ryan Harvey, 2018
            debugging_fig=0;
            
            p = inputParser;
            addRequired(p,'ratemap',@isnumeric)
            addParameter(p,'minPeakRate',2,@isnumeric)
            addParameter(p,'minFieldWidth',2,@isnumeric)
            addParameter(p,'maxFieldWidth',30,@isnumeric)
            addParameter(p,'percentThreshold',.2,@isnumeric)
            
            
            parse(p,varargin{:})
            
            ratemap = p.Results.ratemap;
            minPeakRate = p.Results.minPeakRate;
            minFieldWidth = p.Results.minFieldWidth;
            maxFieldWidth = p.Results.maxFieldWidth;
            threshold = p.Results.percentThreshold; % change between peak rate and start/stop of field
            
            
            warning off  % findpeaks.m throws warnings if peak isn't found...
            
            for i=1:size(ratemap,1)
                fields{i} = [];
                [pks locs w] = findpeaks(fastrms(ratemap(i,:),5),'minpeakheight',minPeakRate,'MinPeakWidth',minFieldWidth);
                exclude=[];
                for j=1:length(locs)-1
                    if min(ratemap(i,locs(j):locs(j+1))) > ((pks(j)+pks(j+1))./2) * threshold
                        % exclude fields without a 90 % decrease in rate between peaks
                        if pks(j) > pks(j+1)
                            exclude = [exclude;j+1];
                        elseif pks(j) < pks(j+1)
                            exclude = [exclude;j];
                        end
                    end
                end
                pks(exclude) = [];
                locs(exclude)=[];
                fieldCount = 1;
                for j=1:length(locs)
                    
                    Map_Field = ratemap(i,:) > pks(j) * threshold;
                    
                    start = locs(j);
                    stop = locs(j);
                    while Map_Field(start) == 1  && start > 1
                        start = start-1;
                    end
                    while Map_Field(stop) == 1  && stop < length(Map_Field) -1
                        stop = stop+1;
                    end
                    if stop - start > minFieldWidth && stop - start < maxFieldWidth
                        fields{i}{fieldCount}.start = start;
                        fields{i}{fieldCount}.stop = stop;
                        fields{i}{fieldCount}.width = stop - start;
                        fields{i}{fieldCount}.peakFR = pks(j);
                        fields{i}{fieldCount}.peakLoc = locs(j);
                        com = start; % calculate center of mass for field
                        fields{i}{fieldCount}.COM = fields{i}{fieldCount}.peakLoc;
                        while sum(ratemap(i,start:stop)) - sum(ratemap(i,start:com)) > sum(ratemap(i,start:stop))./2
                            fields{i}{fieldCount}.COM = com;
                            com = com + 1;
                        end
                        fieldCount = fieldCount + 1;
                    end
                end
                
                % if the peak rate is below 1hz...
                if isempty(fields{i})
                    [fields{i}{1}.peakFR,fields{i}{1}.peakLoc]=max(ratemap(i,:));
                    fields{i}{1}.width=length(ratemap(i,:));
                    fields{i}{1}.start=1;
                    fields{i}{1}.stop=length(ratemap(i,:));
                    fields{i}{1}.COM = fields{i}{1}.peakLoc;
                    com=1;
                    while sum(ratemap(i,fields{i}{1}.start:fields{i}{1}.stop))...
                            - sum(ratemap(i,fields{i}{1}.start:com)) > sum(ratemap(i,fields{i}{1}.start:fields{i}{1}.stop))./2
                        fields{i}{1}.COM = com;
                        com = com + 1;
                    end
                end
                
                % remove fields with the same field boundaries and keep the one with
                % the highest peak rate
                for fie=1:length(fields{i})
                    PR(fie)=fields{i}{fie}.peakFR;
                    start_(fie)=fields{i}{fie}.start;
                    stop_(fie)=fields{i}{fie}.stop;
                end
                fielddets=[start_',stop_',PR'];
                [~,idx]=sort(fielddets(:,3),'descend');
                fielddets=fielddets(idx,:);
                
                [C,ia,~]=unique(fielddets(:,1:2),'rows','stable');
                
                Z=zeros(size(fielddets,1),2);
                Z(ia,:)=C;
                Z=Z(idx,:);
                
                fields_to_delete=find(all(Z==0, 2));
                for f=1:length(fields_to_delete)
                    fields{i}{fields_to_delete(f)}=[];
                end
                fields{i}=fields{i}(~cellfun('isempty',fields{i}));
                
                if debugging_fig
                    for f=1:length(fields{i})
                        figure;
                        plot(ratemap(i,:),'k')
                        grid on
                        hold on
                        plot(fields{i}{f}.start,ratemap(i,fields{i}{f}.start),'*r' )
                        text(fields{i}{f}.start,ratemap(fields{i}{f}.start),'start')
                        
                        plot(fields{i}{f}.stop,ratemap(i,fields{i}{f}.stop),'*r' )
                        text(fields{i}{f}.stop,ratemap(i,fields{i}{f}.stop),'stop')
                        
                        plot(fields{i}{f}.COM,ratemap(i,fields{i}{f}.COM),'*r' )
                        text(fields{i}{f}.COM,ratemap(i,fields{i}{f}.COM),'COM')
                        
                        plot(fields{i}{f}.peakLoc,ratemap(i,fields{i}{f}.peakLoc),'*r' )
                        text(fields{i}{f}.peakLoc,ratemap(i,fields{i}{f}.peakLoc),'peak')
                    end
                end
            end
            warning on
        end
        
        function [ ThPrecess ] = PHprecession(phase,spks_VEL,occ4Ph,fieldbound)
            %PHprecession filters for theta and calculates phase precession
            %
            % Input:    EEG_DownSampledData:        raw downsampled lfp data
            %           EEG_DownSampledTimestamps:  timestamps associated with downsampled lfp data
            %           spks_VEL:                   timestamps associated with spike occurrences
            %           NewSFreq:                   downsampled frequency
            %           track_length:               length or diameter of apparatus
            % ----------------------------------------------------------------------------------------
            % Output:   ThPrecess
            %               phaselock:
            %                           Rlength:R length of spike phases
            %                           Pval:   p value of R length
            %               slope:              slope of regression line
            %               RSquared:           R-Squared
            %               scatteredPH:        phase by position matrix for scatter plot
            %               lapSlope:           mean slope for each lap through firing field
            %               lapR2:              mean R-Squared for each lap
            %               lapCorrelation:     mean correlation for each lap (Kempter et al. 2012)
            %               lapPhaseOffset:     mean phase offset foe each lap (Kempter et al. 2012)
            %               circLinCorr:        Correlation between position and phase (Kempter et al. 2012)
            %               pval:               p value for circ lin correlation (Kempter et al. 2012)
            %               slopeCpU:           slope of regression line in degrees calculated by Kempter method (Kempter et al. 2012)
            %               phaseOffset:        phase offset (Kempter et al. 2012)
            %               RR:                 R-length (Kempter et al. 2012)
            
            % ----------------------------------------------------------------------------------------
            % Ryan E Harvey April 2017;
            % edited April 25th 2018;
            % edited Dec 2nd 2018: to accept multiple fields
            %
            % ADD FMA TO PATH
            
            ThPrecess.phaselock.Rlength=NaN;
            ThPrecess.phaselock.Pval=NaN;
            ThPrecess.slope=NaN;
            ThPrecess.RSquared=NaN;
            ThPrecess.scatteredPH=NaN;
            ThPrecess.lapSlope=NaN;
            ThPrecess.lapR2=NaN;
            ThPrecess.lapCorrelation=NaN;
            ThPrecess.lapPhaseOffset=NaN;
            ThPrecess.circLinCorr=NaN;
            ThPrecess.pval=NaN;
            ThPrecess.slopeCpU=NaN;
            ThPrecess.phaseOffset=NaN;
            ThPrecess.stats=NaN;
            ThPrecess.data=NaN;
            
            if size(spks_VEL,1)<10
                return
            end
            
            % NORM POSITION
            position=[occ4Ph(:,1),rescale(occ4Ph(:,2),0,1),rescale(occ4Ph(:,3),0,1)];
            
            % CHECK TO SEE IF AT LEAST 10 SPIKES ARE WITHIN THE FIELD BOUNDARIES
            [ts,idx]=unique(position(:,1));
            rescalexspk=interp1(ts,position(idx,2),spks_VEL(:,1));
            if sum(rescalexspk>fieldbound(1) & rescalexspk<fieldbound(2))<10
                return
            end
            
            % RUN FMA PRECESSION CODE
            [data,stats]=PhasePrecession(position(idx,:),spks_VEL(:,1),phase,'boundaries',fieldbound);
            ThPrecess.stats=stats;
            ThPrecess.data=data;
            
            
            % PLOT
            % figure;
            % PlotPhasePrecession(data,stats)
            % REMOVE FMA
            
            spks_VEL_working = interp1(phase(:,1),phase(:,2),spks_VEL(:,1)','linear');
            
            % COMPUTE PHASE LOCKING
            ThPrecess.phaselock.Rlength=circ_r(spks_VEL_working');
            [ThPrecess.phaselock.Pval,~]=circ_rtest(spks_VEL_working');
            
            % PLACE OUTPUT IN STRUCTURE
            ThPrecess.slope=stats.slope;
            ThPrecess.RSquared=stats.r2;
            ThPrecess.scatteredPH=[data.position.x,wrapTo360(rad2deg(data.position.phase))];
            lapslope=stats.lap.slope;
            lapslope(isnan(lapslope))=[];
            lapslope(isinf(lapslope))=[];
            ThPrecess.lapSlope=nanmean(lapslope);
            
            lapr2=stats.lap.r2;lapr2(isnan(lapr2))=[];lapr2(isinf(lapr2))=[];
            ThPrecess.lapR2=nanmean(lapr2);
            
            circ_lin_corr=[];
            phi0_deg=[];
            for i=1:length(stats.lap.slope)
                if length(data.position.x(data.position.lap==i))>=2
                    [circ_lin_corr(i,1),pval(i,1),slope_deg(i,1),phi0_deg(i,1)]=...
                        ALLFUNCS.kempter_lincirc(data.position.x(data.position.lap==i),...
                        rad2deg(data.position.phase(data.position.lap==i)));
                else
                    circ_lin_corr(i,1)=NaN;pval(i,1)=NaN;slope_deg(i,1)=NaN;phi0_deg(i,1)=NaN;RR(i,1)=NaN;
                end
            end
            ThPrecess.lapCorrelation=nanmean(circ_lin_corr);
            ThPrecess.lapPhaseOffset=nanmean(phi0_deg);
            
            [ThPrecess.circLinCorr,ThPrecess.pval,ThPrecess.slopeCpU,ThPrecess.phaseOffset]=ALLFUNCS.kempter_lincirc(data.position.x,data.position.phase);
        end
        
        function [ rho,p,s,b ] = kempter_lincirc( x,theta,varargin )
            %KEMPTER_LINCIRC - Linear-circular correlation
            %   Performes a linear-circular correlation (Kempter et al, 2012)
            %
            %   [ RHO,P,S,B ] = kempter_lincirc(X,THETA)
            %   [ RHO,P,S,B ] = kempter_lincirc(X,THETA,S,B)
            %
            %   ARGUMENTS
            %   * X: A vector of linear values
            %   * THETA: A vector, the same size as X, containing circular values (in
            %   radians)
            %   * S: The presumed best slope in cycles/unit. If not entered, uses
            %   anglereg
            %   * B: The presumed best phase shift in radians. If not entered, uses
            %   anglereg
            %
            %   RETURNS
            %   * rho: The linear-circular correlation coefficient
            %   * p: The statical significance of the correlation
            %   * S: The slope in cycles per unit. To convert to radians, multiply by
            %   2*pi
            %   * B: The phase offset in radians.
            %
            % This code has been freely distributed by the authors. If used or
            % modified, we would appreciate it if you cited our paper:
            % Climer, J. R., Newman, E. L. and Hasselmo, M. E. (2013), Phase coding by
            %   grid cells in unconstrained environments: two-dimensional phase
            %   precession. European Journal of Neuroscience, 38: 2526?2541. doi:
            %   10.1111/ejn.12256
            %
            % RELEASE NOTES
            %   v1.0 2014-10-15 Release (Jason Climer, jason.r.climer@gmail.com)
            %
            % This file is part of pass_index. All or part of this file may be
            % considered derivative.
            %
            % Copyright (c) 2014, Trustees of Boston University
            % All rights reserved.
            %
            % Redistribution and use in source and binary forms, with or without
            % modification, are permitted provided that the following conditions are
            % met:
            %
            % 1. Redistributions of source code must retain the above copyright notice,
            % this list of conditions and the following disclaimer.
            %
            % 2. Redistributions in binary form must reproduce the above copyright
            % notice, this list of conditions and the following disclaimer in the
            % documentation and/or other materials provided with the distribution.
            %
            % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
            % "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
            % TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
            % PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
            % OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
            % EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
            % PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
            % PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
            % LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
            % NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
            % SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
            
            % Parse inputs
            if nargin==1, s=varargin{1}; varargin = {}; end;
            if strcmp(class(x),'single')
                x = double(x);
            end
            
            goods = ~isnan(x)&~isnan(theta);
            x = x(goods);
            theta = theta(goods);
            
            % Angular regression
            if ~exist('s','var')
                [s,b]=ALLFUNCS.anglereg(x,theta,varargin{:});
            end
            
            n = length(x);
            
            % Calculate rho
            phi = mod(s*x,2*pi);
            theta = mod(theta,2*pi);
            phi_ = angle(sum(exp(1i*phi))/n);
            theta_ = angle(sum(exp(1i*theta))/n);
            rho = abs(sum(sin(theta-theta_).*sin(phi-phi_))/...
                sqrt(sum(sin(theta-theta_).^2)*sum(sin(phi-phi_).^2)))*sign(s);
            
            % Calculate p
            lambda = @(i,j)n^-1*sum((sin(phi-phi_).^i).*(sin(theta-theta_).^j));
            z = rho*sqrt(n*lambda(2,0)*lambda(0,2)/lambda(2,2));
            p = 1-erf(abs(z)/sqrt(2));
        end
        
        function [ s,b ] = anglereg( x, theta, bnds)
            % ANGLEREG - Linear-circular regression
            %   Performes a linear-circular regression (Kempter et al, 2012)
            %
            %   [ S,B ] = ANGLEREG(X,THETA)
            %   [ S,B ] = ANGLEREG(X,THETA,TOLERANCE)
            %   [ S,B ] = ANGLEREG(X,THETA,TOLERANCE,MAXITER)
            %
            %   ARGUMENTS
            %   * x: A vector of linear values
            %   * theta: A vector, the same size as X, containing circular values (in
            %   radians)
            %   * bnds (optional): A 2 element vector with the lower and upper bounds
            %   of the slope (B)
            %
            %   RETURNS
            %   * S: The slope in cycles per unit. To convert to radians, multiply by
            %   2*pi
            %   * B: The phase offset in radians.
            %
            % This code has been freely distributed by the authors. If used or
            % modified, we would appreciate it if you cited our paper:
            % Climer, J. R., Newman, E. L. and Hasselmo, M. E. (2013), Phase coding by
            %   grid cells in unconstrained environments: two-dimensional phase
            %   precession. European Journal of Neuroscience, 38: 2526?2541. doi:
            %   10.1111/ejn.12256
            %
            % RELEASE NOTES
            %   v1.0 2014-10-15 Release (Jason Climer, jason.r.climer@gmail.com)
            %
            % This file is part of pass_index. Parts of this file may be
            % considered derivative.
            %
            % Copyright (c) 2014, Trustees of Boston University
            % All rights reserved.
            %
            % Redistribution and use in source and binary forms, with or without
            % modification, are permitted provided that the following conditions are
            % met:
            %
            % 1. Redistributions of source code must retain the above copyright notice,
            % this list of conditions and the following disclaimer.
            %
            % 2. Redistributions in binary form must reproduce the above copyright
            % notice, this list of conditions and the following disclaimer in the
            % documentation and/or other materials provided with the distribution.
            %
            % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
            % "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
            % TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
            % PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
            % OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
            % EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
            % PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
            % PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
            % LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
            % NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
            % SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
            
            % Parse inputs
            if ~exist('bnds','var')
                bnds = [];
            end
            
            % format inputs
            theta = mod(theta,2*pi);
            if size(x,1)>size(x,2)
                x = x';
            end
            if size(theta,1)>size(theta,2)
                theta = theta';
            end
            
            X = [ones(size(x(:))) x(:)];
            
            
            
            if isempty(bnds)% No bounds
                % Initial guess from least squares and rolling phase
                phi = fminbnd(@(phi)sum((mod(theta(:)+phi,2*pi)-X*(inv(X'*X)*X'*mod(theta(:)+phi,2*pi))).^2),0,2*pi,...
                    optimset('Display','off','TolX',1e-10));
                b = inv(X'*X)*X'*mod(theta(:)+phi,2*pi);
                b(1) = b(1)-phi;
                n = length(x);
                
                % Use above initial guess...
                [s1,r1] = fminsearch(...
                    @(s)-sqrt((1/n*sum(cos(theta-2*pi*s*x)))^2+(1/n*sum(sin(theta-2*pi*s*x)))^2),...
                    b(2));
                % And the perpendicular slope....
                [s2,r2] = fminsearch(...
                    @(s)-sqrt((1/n*sum(cos(theta-2*pi*s*x)))^2+(1/n*sum(sin(theta-2*pi*s*x)))^2),...
                    -1/b(2)...
                    ,optimset('Display','off','TolX',1e-100,'TolFun',1e-100));
                % And take the better one
                if -r1>-r2
                    s = s1;
                else
                    s = s2;
                end
            else
                % Use the bounds the user has given
                s = fminbnd(@(s)-sqrt((1/n*sum(cos(theta-2*pi*s*x)))^2+(1/n*sum(sin(theta-2*pi*s*x)))^2),bnds(1),bnds(2));
            end
            
            % Slope
            b = atan2(sum(sin(theta-2*pi*s*x)),sum(cos(theta-2*pi*s*x)));
        end
        
        function [IC]=infocontent(SmoothRateMap,occ)
            for i=1:size(SmoothRateMap,1)
                [nBinsy,nBinsx]=size(occ(i,:));
                
                rY=reshape(SmoothRateMap(i,:),nBinsx*nBinsy,1);
                rY(isnan(rY))=0;
                rY(isinf(rY))=0;
                occRSHP=reshape(occ(i,:),nBinsx*nBinsy,1);
                occSUM=sum(occRSHP);
                pX=occRSHP./occSUM;
                [nBins,nCells]=size(rY);
                relR=rY./kron(ones(nBins,1),pX'*rY);
                log_relR=log2(relR);
                log_relR(isinf(log_relR))=0;
                IC(i)=sum(kron(pX,ones(1,nCells)).*relR.*log_relR);
            end
        end
        
        function [trialmaps,occ,frames_w_spk,trial_frames_w_spk]=split_trials(data,spikes)
            frames_w_spk=[];
            for i=1:size(data.trial_ts,1)
                tempframes=data.frames(data.frames(:,1)>data.trial_ts(i,1) & data.frames(:,1)<data.trial_ts(i,2),:);
                
                tempspikes=spikes(spikes>data.trial_ts(i,1) & spikes<data.trial_ts(i,2));
                [spkmatrix]=ALLFUNCS.get_spkmatrix(tempframes,tempspikes);
                
                trial_frames_w_spk{i}=sortrows([[spkmatrix,ones(size(spkmatrix,1),1)];...
                    [tempframes,zeros(length(tempframes),1)]],1);
                
                % get frames with embedded spike positions & spike binary
                frames_w_spk=[frames_w_spk;sortrows([[spkmatrix,ones(size(spkmatrix,1),1)];...
                    [tempframes,zeros(length(tempframes),1)]],1)];
                
                [trialmaps(i,:),occ(i,:)]=ALLFUNCS.binlaps(data.frames,tempframes,data.samplerate,spkmatrix,data.mazesize);
            end
            trialmaps(isnan(trialmaps))=0;
            trialmaps(occ==0)=NaN;
            
            frames_w_spk=sortrows(frames_w_spk,1);
        end
        
        function [right,left]=split_trials_RL(ratemaps,occ,trial_frames_w_spk)
            right_idx=1:2:size(ratemaps,1);
            left_idx=2:2:size(ratemaps,1);

            right.ratemaps=ratemaps(right_idx,:);
            right.occ=occ(right_idx,:);
            right.frames_w_spk=cat(1,trial_frames_w_spk{right_idx});
            right.trial_frames_w_spk={trial_frames_w_spk{right_idx}};

            
            left.ratemaps=ratemaps(left_idx,:);
            left.occ=occ(left_idx,:);
            left.frames_w_spk=cat(1,trial_frames_w_spk{left_idx});
            left.trial_frames_w_spk={trial_frames_w_spk{left_idx}};
        end
    end
end