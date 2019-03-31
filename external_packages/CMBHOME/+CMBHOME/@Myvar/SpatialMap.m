function [rate_map, xdim, ydim, occupancy, no_occupancy] = SpatialMap(self,variableName,root,varargin)
% [rate_map, xs, ys] = root.plot_rate_map(cel, varargin)
%
% Plots spatial rate map for cell cel = [ tetrode, cell ]
%
% Returns matrix rate_map with x and y dimensions in vectors xs and ys. 
%
% ARGUMENTS
%   cel             1 x 2 vector like [tetrode index, cell index]
%   params          (see below)
%
% RETURNS
%   rate_map        matrix of F/bin (if continuize_epochs == 0), a third
%                   dimension exists for each epoch
%   xdim            vector of x dimensions
%   ydim            vector of y dimensions
%   occupancy       matrix the size of rate_map indicating occupancy in
%                   seconds
%   no_occupancy    logical matrix the size of rate_map. 1's indicate no
%                   occupancy. useful for AlphaData property of imagesc
%
% OPTIONAL PARAMETERS
%   
%   xdim                vector of bin edges along x dimension
%   ydim                vector of bin edges along y dimension
%   continuize_epochs   0 or 1 (0). If 0, ratemap is calculated for each
%                       epoch (adds 3rd dim to rate_map output, if 1, ratemap
%                       is calculated across all epochs
%   supress_plot        0 or 1 (1). If 0, plots ratemap
%   figure_handle       If supplied, plot the rate map to the given figure
%   std_smooth_kernel   STD (cm) of the gaussian kernel to smooth the rate map  
%   binside             The length in cm of the side of a bin when
%                       calculating the rate map
%   show_scale          (1) if plotting, shows scale or not
%   show_peak_f         (1) if 1, shows peak f or not
%   show_mean_f         (1) if 1, shows mean f or not
%   omit_islands        (1) if 1, solves for center of mass of occupancy,
%                       and then looks for any stray pixels and sets their 
%                       spikes to zero  
%   omit_noocc          (1) if 1, sets to zero pixels from ratemap which had
%                       occupancy values = 0

% wchapman 20130730

import CMBHOME.Utils.*

p = inputParser;

p.addParamValue('xdim', [], @isnumeric);
p.addParamValue('ydim', [], @isnumeric);
p.addParamValue('clims', [], @(c) numel(c)==2);
p.addParamValue('continuize_epochs', 0, @(c) numel(c)==1 && (c==1 || c==0));
p.addParamValue('ifplot', 0, @(c) numel(c)==1 && (c==1 || c==0));
p.addParamValue('figure_handle', '', @(c) numel(c)==1);
p.addParamValue('std_smooth_kernel', 3, @isnumeric);
p.addParamValue('binside', 3, @isnumeric)
p.addParamValue('show_scale', 1, @isnumeric)
p.addParamValue('show_peak_f', 1, @isnumeric)
p.addParamValue('show_mean_f', 1, @isnumeric)
p.addParamValue('omit_islands', 1, @isnumeric)
p.addParamValue('omit_noocc', 1, @isnumeric)

p.parse(varargin{:});

xdim = p.Results.xdim;
ydim = p.Results.ydim;
clims = p.Results.clims;
continuize_epochs = p.Results.continuize_epochs;
ifplot = p.Results.ifplot;
figure_handle = p.Results.figure_handle;
std_smooth_kernel = p.Results.std_smooth_kernel;
binside = p.Results.binside;
show_scale = p.Results.show_scale;
show_peak_f = p.Results.show_peak_f;
show_mean_f = p.Results.show_mean_f;
omit_islands = p.Results.omit_islands;
omit_noocc = p.Results.omit_noocc;


%if isnan(cel), cel = root.cel; end

%if any(size(cel)~=[1, 2]), error('cell must be size 1x2 like [tetrode, cell]'); end

if continuize_epochs, root = MergeEpochs(root); end % so that we don't double count data

[occupancy, xdim, ydim] = root.Occupancy(xdim, ydim, continuize_epochs, binside);
occupancy= occupancy*root.fs_video; %bring back to samples
binside = root.spatial_scale * diff(xdim(1:2)); % solve for binside in case xdim and ydim were specified

if omit_islands, occupancy = OmitIslands(occupancy); end

no_occupancy = occupancy==0; % mark indeces where there was no occupancy so we can correct after smoothing
    
if (continuize_epochs && iscell(root.ind))
    [var] = ContinuizeEpochs(self.get(variableName)); 
else
    var = self.get(variableName);
end

isCirc = self.isCirc_get(variableName);

if ~iscell(var)
    %spikes = hist3([spk_x, spk_y], 'Edges', {xdim, ydim});
    [~,wx] = histc(root.x,xdim);
    [~,wy] = histc(root.y,ydim);
    varSum = NaN(length(xdim),length(ydim));
    
    if ~isCirc
        for x = 1:length(xdim)
            for y = 1:length(ydim)
                inds = (wx==x) & (wy==y);
                varSum(x,y) = nansum(var(inds));
            end
        end
        rate_map = SmoothMat(varSum, [5*std_smooth_kernel/binside, 5*std_smooth_kernel/binside], std_smooth_kernel/binside) ./SmoothMat(occupancy, [5*std_smooth_kernel/binside, 5*std_smooth_kernel/binside], std_smooth_kernel/binside); % smooth the spikes and occupancy with a 5x5 bin gaussian with std=1

    else
        for x = 1:length(xdim)
            for y = 1:length(ydim)
                inds = (wx==x) & (wy==y);
                varSum(x,y) = nansum(exp(1i*var(inds)));
            end
        end
        rate_map = SmoothMat(varSum, [5*std_smooth_kernel/binside, 5*std_smooth_kernel/binside], std_smooth_kernel/binside) ./SmoothMat(occupancy, [5*std_smooth_kernel/binside, 5*std_smooth_kernel/binside], std_smooth_kernel/binside); % smooth the spikes and occupancy with a 5x5 bin gaussian with std=1
        rate_map = angle(rate_map);
    end
        
    if omit_noocc
        rate_map(no_occupancy) = 0; % set no occupancy to zero
    end

    rate_map = rate_map'; % returns these three
    
    occupancy = occupancy';

    no_occupancy = no_occupancy';
    
else % multiple epochs

    rate_map = zeros(size(occupancy, 2), size(occupancy, 1), size(occupancy, 3));

    new_occupancy = zeros(size(occupancy, 2), size(occupancy,1), size(occupancy,3));

    for i = 1:size(occupancy,3)
        vt = var{i};
        x = root.x{i}; y = root.y{i};
        [~,wx] = histc(x,xdim);
        [~,wy] = histc(y,ydim);
        varSum = NaN(length(xdim),length(ydim));
        for x = 1:length(xdim)
            for y = 1:length(ydim)
                inds = (wx==x) & (wy==y);
                varSum(x,y) = nansum(vt(inds));
            end
        end

        tmp = SmoothMat(varSum, [5*std_smooth_kernel/binside, 5*std_smooth_kernel/binside], std_smooth_kernel/binside)./SmoothMat(occupancy(:,:,i), [5*std_smooth_kernel/binside, 5*std_smooth_kernel/binside], std_smooth_kernel/binside);

        if omit_noocc
            tmp(no_occupancy(:,:,i)') = 0;
        end

        rate_map(:, :, i) = tmp';

        new_occupancy(:,:,i) = no_occupancy(:,:,i)';

    end

    occupancy = permute(occupancy,[2 1 3]);
    
    no_occupancy = new_occupancy;

end
        
if ifplot && (size(rate_map,3)==1 || ndims(rate_map)==2)
    PlotIt(rate_map, no_occupancy, clims, xdim, ydim, figure_handle, root, show_scale, show_peak_f, show_mean_f, nanmean(rate_map(:)));
end
    
end

function occupancy = OmitIslands(occupancy)
% Takes matrix of occupancy values, calculates center of mass of pixels>0,
% and then finds all disconnected pixels and sets them = 0

%subplot(1, 2, 1), imagesc(occupancy)

s = regionprops(occupancy>0, {'FilledArea', 'PixelIdxList'});

l = numel(s);

areas = vertcat(s.FilledArea);

[~, inds] = sort(areas);

for i = 1:length(inds)-1
    
    occupancy(s(inds(i)).PixelIdxList) = 0;
    
end 

end

function PlotIt(rate_map, no_occupancy, clims, xdim, ydim, handle, root, show_scale, show_peak_f, show_mean_f, mean_f)

    if all(isnan(rate_map)) | all(rate_map==0)
        
        text(.05, .4, 'No figure'), axis off
        
        return
        
    end

    xs = [min(xdim) max(xdim)];
    ys = [min(ydim) max(ydim)];

    if ~isempty(handle), figure(handle); end % make current the figure passed

    if isempty(clims)
        clims = [min(min(rate_map)) max(max(rate_map))];
    end
    
    pad = [-.03 .02]; % percent pad plot
    
    [cbar, clims] = CMBHOME.Utils.SmartColorbar(clims, 'jet(255)');

    rate_map(no_occupancy) = clims(1);

    imagesc(xdim, ydim, rate_map, clims); hold on;
    
    colormap(cbar);
    
    axis equal off

    xlim(diff(xs).*pad+xs);
    ylim(diff(ys).*pad+ys);
    
    if show_scale
        line([xs(1)+.75*diff(xs), xs(2)], [-.03*diff(ys)+ys(1), -.03*diff(ys)+ys(1)], 'Color', 'k', 'LineWidth', 2);
        text(xs(2), -.03*diff(ys)+ys(1), [num2str(.25*diff(xs)*root.spatial_scale, 3) ' cm'], 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlign', 'right', 'VerticalAlign', 'bottom');
    end
    
    str_f = '';
    
    if show_peak_f, str_f = ['p: ' num2str(max(max(rate_map)), 2) 'Hz']; end
    
    if show_mean_f, str_f = [str_f ' m: ' num2str(mean_f, 2) 'Hz']; end
    
    if show_peak_f | show_mean_f
        text(xs(2), .045*diff(ys)+ys(2), str_f, 'FontSize', 6.8, 'FontWeight', 'bold', 'HorizontalAlign', 'right');
    end
        
    rmpos=get(gca, 'Position');

    set(gca,'YDir','normal'); % so plotting functions dont reverse axis

end
