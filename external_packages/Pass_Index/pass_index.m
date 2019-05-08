function results = pass_index( varargin )
% PASS_INDEX - Calcultes the pass index and plots
%
% Calculates the pass index for the data passed. This function uses the
% pass_index_parser to generate its input structure.
%
% RESULTS = PASS_INDEX_PARSER(POS_TS,POS,SPK_TS,LFP_TS,LFP_SIG)
% RESULTS = PASS_INDEX_PARSER(POS_TS,POS,SPK_TS,LFP_TS,LFP_SIG,PARAMS)
%
%   ARGUMENTS
%   * POS_TS: Vector of time stamps for the sample state
%   * POS: MXN matrix of the sample state, where M is the number of samples
%   and N is the dimensions of POS
%   * SPK_TS: Spike times for the cell
%   * LFP_TS: Time stamps for the local field potential (LFP) Sample
%   * LFP_SIG: The LFP signal
%
%   PARAMETERS
%   * plots: Default false. If false, doesn't plot. If true or 'all', plots
%   all possible plots. Can be true,'all',any or or a cell array of the
%   following:
%     	Trajectory: (1 to 3D only) This plots the trajectory of the animal with the spikes of the cell colored by the pass index of the spike.
%       Rate map: (1 to 3D only) This plots the rate map of the cell.
%   	Field index map: (1 to 3D only) This plots the field index map of the cell. See table 1 for notes on custom implementation.
%   	Scatter plot: Shows the pass index versus two cycles of the lfp phase for all spikes, and calculates the linear circular correlation using the circular-linear correlation19. This can be calculated using the [corr_val,p,s,b]=kempter_lincirc(linear,circular) function included, and are output in the results struct (see 5.4 and Table 1).
%   	Density map: Shows the density of the scatter plot, using 100 phase bins and 40 pass index bins, then smoothing with a pseudo-Gaussian kernel with width of 1.5 pixel.
%   * subplots: Default []. If the parameter ?subplots? is set to a nX2
%   matrix, where n is the number of plots, each of the rows of the
%   subplots parameter will  be used to determine where the subplot will
%   go. Otherwise, will plot as squarly as possible, biased wider than
%   taller.
%
%   PASS INDEX PARAMETERS
%   * method: Default 'grid'. Can be 'grid','place', or custom. Updates
%   other unset fields for these techniques.
%   * binside: Default 2*N, where N is the dimensionality of POS. Side of
%   the bins for rate mapping.
%   * smth_width: Default 3*BINSIDE, width of Gaussian smoothing kernel
%   * field_index: Default @field_index_fun, can be a vector of the same
%   number of elements as pos_ts, or can be a function handle which takes
%   in the same parameters as pass_index.
%   * sample_along: Default 'auto', can be 'arc_length', 'raw_ts', or a
%   nX2 matrix where n is the number of resampled steps, the first column
%   is the resampled timestamps and the second column is the sampled values,
%   or a function handle that returns a nX2 matrix as described above. Set
%   from 'auto' to 'arc_length' if method is 'place' or 'grid'.
%   * filter_band: Default 'auto', can be any positive frequency range in
%   cycles/unit sampled along using the ?filter_band? parameter.
%   Additionally, filter_band can be a function handle which returns a
%   modified signal. Set from 'auto' to [0.0749 0.0029] if 'method' is
%   'grid' and to the [3*D 1/6*D].^-1, where D is the field width
%   determined by finding the N-dimensional volume of the region with at
%   least 10% of the maximum firing rate, and calculating the diameter of
%   the n-ball with the same volume.
%   * lfp_filter: Default [6 10]. can be changed to any frequency range in
%   Hz as [low high] or as a function handle with the form lfp_phases =
%   custom_phase_func(lfp_ts,lfp_sig) for custom phase estimation, for
%   example, by taking asymmetry into account
%
%   RETURNS
%   * RESULTS: A struct containing the following fields:
%       rate_map	Occupancy normalized rate map for the cell.
%       field_index_map	Field index map as calculated by the field index function.
%       centers the centers for the maps.
%       ts	Resampled time stamps.
%       cs	Resampled axis positions. As a default, the distance along the trajectory (cm).
%       field_index	Field index at every point along the resampled axis.
%       filtered_field_index Field index after filter is applied
%       pass_index	Pass index at every point along the resampled axis.
%       spk_pass_index	The pass index at the time of each spike.
%       spk_theta_phase	The theta phase at the time of each spike.
%       rho	The linear-circular correlation coefficient from kempter_lincirc.
%       p	The significance level of the correlation.
%       s	The slope in ? per pass.
%       is_precessing	True if p<0.05 and -1440<s<-22 ? per pass.
%       density	Density map, with LFP phase broken into 100 bins and pass index broken into 40.
%       filtered_lfp	The filtered LFP
%       filtered_lfp_phase	The phase of the filtered LFP
%
% This code has been freely distributed by the authors. If used or
% modified, we would appreciate it if you cited our paper:
% Climer, J. R., Newman, E. L. and Hasselmo, M. E. (2013), Phase coding by 
%   grid cells in unconstrained environments: two-dimensional phase 
%   precession. European Journal of Neuroscience, 38: 2526–2541. doi: 
%   10.1111/ejn.12256
%
% RELEASE NOTES
%   v1.0 2014-10-15 Release (Jason Climer, jason.r.climer@gmail.com)
%
% This file is part of pass_index.
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
P = pass_index_parser(varargin{:});
for i = fields(P.Results)'
    eval([i{1} ' = P.Results.' i{1} ';']);
end

% Calculate field index
if isequal(class(field_index),'function_handle')
    fi_fun = field_index;
    field_index = fi_fun(varargin{:});
    if any(cellfun(@(x)isequal(x,'field_index'),varargin))
        varargin{find(cellfun(@(x)isequal(x,'field_index'),varargin))+1}=field_index;
    else
        varargin = [varargin {'field_index',field_index}];
    end
end

% Resample
if isequal(class(sample_along),'function_handle')
    sample_along = sample_along(varargin{:});
    if any(cellfun(@(x)isequal(x,'sample_along'),varargin))
        varargin{find(cellfun(@(x)isequal(x,'sample_along'),varargin))+1}=sample_along;
    else
        varargin = [varargin {'sample_along',sample_along}];
    end
end

% Filter
[filtered] = filter_band(varargin{:});

% Calculate Pass Index
pass_index = angle(hilbert(filtered))/pi;
spk_pass_index = (mod(interp1(sample_along(:,2),unwrap(pass_index*pi),spk_ts,'nearest','extrap')+pi,2*pi)-pi)/pi;

% Calculate LFP
if ~isempty(lfp_sig)
    [filtered_lfp,lfp_phase] = lfp_filter(varargin{:});
    spk_theta_phase = mod(interp1(lfp_ts,unwrap(lfp_phase),spk_ts)+pi,2*pi)-pi;
end

%% Packing results
results = struct();

[map,centers,occupancy] = rate_map(varargin{:});
results.rate_map = map;
results.centers = centers;
results.occupancy = occupancy;

fi_map = NaN;

try
    fi_map = fi_fun(varargin{:},'get_map',true);
    if ~isequal(size(fi_map),size(map))
        fi_map = NaN;
    end
catch err
end

results.field_index_map = fi_map;
results.ts = sample_along(:,2);
results.cs = sample_along(:,1);
results.field_index = sample_along(:,3);
results.pass_index = pass_index;
results.filtered_field_index = filtered;
results.spk_pass_index = spk_pass_index;
results.spk_theta_phase = spk_theta_phase;
results.filtered_lfp = filtered_lfp;
results.filtered_lfp_phase = lfp_phase;

[rho,p,s,b] = kempter_lincirc(spk_pass_index,spk_theta_phase); % Correlation
results.rho = rho;
results.p = p;
results.s = s;
results.b = b;
results.is_precessing = p<0.05&&rad2deg(pi*s)<-22&&rad2deg(pi*2)>-1440;

% Density
try
   dens_oc = histcn([interp1(pos_ts,pass_index,lfp_ts,'nearest')' mod(lfp_phase,2*pi)'],linspace(-1,1,41),linspace(0,2*pi,101));
catch err % Needs interpolation before we can calculate occupancy
   % Interpolate & calculate occupancy for the density map
   dens_oc = zeros(41,101);
   interpd = mod(interp1(pos_ts,unwrap(pass_index*pi)/pi,lfp_ts,'nearest')+1,2)-1;% Interpolate unwrapped positions
   [~,pihist] = histc(interpd,linspace(-1,1,41));
   for i=1:41
      dens_oc(i,:) = histc(mod(lfp_phase(pihist==i),2*pi),linspace(0,2*pi,101));
   end
   dens_oc = dens_oc(1:40,1:100);
end

dens_oc = dens_oc*mean(diff(lfp_ts));
% Calculate density map (spike counts)
[density,~,dens_centers] = histcn([spk_pass_index mod(spk_theta_phase,2*pi)],linspace(-1,1,41),linspace(0,2*pi,101));

% ryan H 2019 bug fix: density had an extra row and column of zeros
if size(density)~=size(dens_oc)
    density = density(1:40,1:100);
end

density = density./dens_oc;% Normalize by occupancy
% Smooth kernel
h = 1.5;
myfilter = fspecial('gaussian',[4 4]*h, h);
density = imfilter(density,myfilter,'replicate');
results.density = density;

%% Plotting
% Parse plot stuff
all_plots = {'Trajectory','Rate map','Field index map','Scatter plot','Density map'};
ip = P;
ip.KeepUnmatched=true;
ip.addParamValue('plots',{},@(x)(iscellstr(x)&&...
    all(ismember(x,all_plots)))||isequal(x,'all')||isequal(x,0)||isequal(x,1)||...
    (ischar(x)&&ismember(x,plots)));
ip.parse(varargin{:});

plots = ip.Results.plots;
if ~iscell(plots)% Not specific plot list
    switch plots
        case {1,'all'}
            plots = all_plots;
        case 0
            plots = {};
        otherwise
            plots = {plots};% Or it is just one plot
    end
end

% Subplots
subplots = floor(sqrt(numel(plots)));
subplots = [subplots ceil(numel(plots)/subplots)];
ip.addParamValue('subplots',subplots,@(x)isequal(x,floor(x))&&numel(x)==2);
ip.addParamValue('units','cm',@ischar);
ip.parse(varargin{:});
subplots = ip.Results.subplots;
units = ip.Results.units;

%% Plot
if numel(plots)>0
    set(gcf,'color','w');
    for i=1:numel(plots)
        subplot(subplots(1),subplots(2),i);
        switch plots{i}
            case 'Trajectory'
                switch size(pos,2)
                    case 1 % 1D
                        plot(pos_ts,pos);
                        xlim(minmax(pos_ts(:)'));
                        ylim([-1 1]*range(pos(:))*1.1+nanmean(pos(:)));
                        
                        xlabel('Time (s)');
                        ylabel(['State (' units ')']);
                        
                        spkpos = spk_pos(pos_ts,pos,spk_ts);
                        hold on;
                        scatter(spk_ts,spkpos,'o','filled','CData',spk_pass_index,'SizeData',20);
                        hold off
                    case 2 %2D
                        plot(pos(:,1),pos(:,2),'k');
                        spkpos = spk_pos(pos_ts,pos,spk_ts);
                        hold on;
                        scatter(spkpos(:,1),spkpos(:,2),'o','filled','CData',spk_pass_index,'SizeData',20);
                        line(...
                            [-floor(0.2*range(pos(:,1))/5)*5 0]+0.95*max(pos(:,1)),...
                            [1 1]*(min(pos(:,2)-0.1*range(pos(:,2)))),...
                            'Color','k','LineWidth',3);
                        text(0.95*max(pos(:,1)),...
                            min(pos(:,2))-0.15*range(pos(:,2)),...
                            [num2str(floor(0.2*range(pos(:,1))/5)*5) ' ' units],...
                            'VerticalAlignment','Cap','HorizontalAlignment','right');
                        axis off equal;
                        hold off;
                        
                    case 3 % 3D
                        plot3(pos(:,1),pos(:,2),pos(:,3),'k','LineWidth',1)
                        spkpos = spk_pos(pos_ts,pos,spk_ts);
                        hold on;
                        scatter3(spkpos(:,1),spkpos(:,2),spkpos(:,3),...
                            'o','filled','CData',spk_pass_index);
                        l = floor(max(range(pos))*0.2/10)*10;
                        line([0 0 0;l 0 0]+min(pos(:,1))-0.01*range(pos(:,1)),...
                            [0 0 0;0 l 0]+min(pos(:,2))-0.01*range(pos(:,2)),...
                            [0 0 0;0 0 l]+min(pos(:,3))-0.01*range(pos(:,3)),...
                            'Color','k','LineWidth',3);
                        
                        text(min(pos(:,1))-0.02*range(pos(:,1)),...
                            min(pos(:,2))-0.02*range(pos(:,2)),...
                            min(pos(:,3))-0.02*range(pos(:,3)),...
                            [num2str(l) ' ' units],...
                            'VerticalAlignment','Cap','HorizontalAlignment','right');
                        hold off;
                        axis off equal;
                        view(3);
                    otherwise % Can't plot
                        warning('pass_index:bad_dim_subplot','Cannot plot over 3d');
                end
            case 'Rate map'
                switch size(pos,2)
                    case 1 % 1D
                        bar(centers{:},map);
                        xlim(minmax(centers{:}));
                        ylim([0 max(map(:))*1.1]);
                        
                        temp = occupancy==0;
                        hold on;
                        while any(temp)
                            j = find(temp,1);
                            k = find(temp(j:end)==0,1)-2;
                            rectangle('Position',[centers{1}(j) 0 range(centers{1}(j:(j+k))) max(map(:))*1.1],'FaceColor',[1 1 1]*0.7,'LineStyle','none');
                            temp(j:(j+k))=0;
                        end
                        hold off;
                        ylabel(['State (' units ')']);
                        xlabel(['Rate (Hz)']);
                    case 2 % 2D
                        pk = quantile(map(occupancy>0),0.99);
                        [cbar, clims] = smart_colorbar([0 pk], jet(255));
                        map(occupancy==0)=clims(1);
                        imagesc(centers{:},map');
                        colormap(cbar);
                        set(gca,'CLim',clims);
                        
                        hold on;
                        line(...
                            [-floor(0.2*range(pos(:,1))/5)*5 0]+0.95*max(pos(:,1)),...
                            [1 1]*(min(pos(:,2)-0.1*range(pos(:,2)))),...
                            'Color','k','LineWidth',3);
                        text(0.95*max(pos(:,1)),...
                            min(pos(:,2))-0.15*range(pos(:,2)),...
                            [num2str(floor(0.2*range(pos(:,1))/5)*5) ' ' units],...
                            'VerticalAlignment','cap','HorizontalAlignment','right');
                        text(max(pos(:,1)),...
                            max(pos(:,2)),...
                            [sprintf('%3.1f',pk) ' Hz'],...
                            'VerticalAlignment','bottom','HorizontalAlignment','right');
                        axis off equal;
                        set(gca,'YDir','normal');
                        hold off;
                        freezeColors;
                    case 3 % 3D
                        set(gcf,'Renderer','OpenGL');
                        pk = quantile(map(occupancy>0),0.99);
                        [cbar, clims] = smart_colorbar([0 pk], jet(255));
                        map(occupancy==0)=clims(1);
                        j = {0 1 2};
                        k = find(occupancy>0)';
                        view(3);
                        for k=k;
                            [j{1} j{2} j{3}] = ind2sub(size(map),k);
                            clr = round((map(k)-min(clims))/range(clims)*size(cbar,1))+1;
                            if clr>size(cbar,1), clr=size(cbar,1); end;
                            plotcube([1 1 1]*binside,...
                                cellfun(@(x,y)x(y),centers,j)-[0.5 0.5 0.5]*binside,...
                                cbar(clr,:),...
                                min((map(k)-min(clims))/range(clims),1)*0.5+0.02,...
                                'EdgeAlpha',0,...
                                'BackfaceCull',1);
                        end
                        l = floor(max(range(pos))*0.2/10)*10;
                        line([0 0 0;l 0 0]+min(pos(:,1))-0.01*range(pos(:,1)),...
                            [0 0 0;0 l 0]+min(pos(:,2))-0.01*range(pos(:,2)),...
                            [0 0 0;0 0 l]+min(pos(:,3))-0.01*range(pos(:,3)),...
                            'Color','k','LineWidth',3);
                        
                        text(min(pos(:,1))-0.02*range(pos(:,1)),...
                            min(pos(:,2))-0.02*range(pos(:,2)),...
                            min(pos(:,3))-0.02*range(pos(:,3)),...
                            [num2str(l) ' ' units],...
                            'VerticalAlignment','Cap','HorizontalAlignment','right');
                        
                        text(min(pos(:,1))-0.02*range(pos(:,1)),...
                            min(pos(:,2))-0.02*range(pos(:,2)),...
                            min(pos(:,3))-0.12*range(pos(:,3)),...
                            [sprintf('%3.1f',pk) ' Hz'],...
                            'VerticalAlignment','bottom','HorizontalAlignment','right');
                        hold off;
                        axis off equal;
                        map(occupancy==0)=clims(1);
                    otherwise % Can't plot
                        warning('pass_index:bad_dim_subplot','Cannot plot over 3d');
                end
            case 'Field index map'
                if ~all(isnan(fi_map(:)))
                    switch size(pos,2)
                        case 1 %1D
                            bar(centers{:},fi_map);
                            xlim(minmax(centers{:}));
                            ylim([0 max(fi_map(:))*1.1]);
                            
                            temp = occupancy==0;
                            hold on;
                            while any(temp)
                                j = find(temp,1);
                                k = find(temp(j:end)==0,1)-2;
                                rectangle('Position',[centers{1}(j) 0 range(centers{1}(j:(j+k))) max(map(:))*1.1],'FaceColor',[1 1 1]*0.7,'LineStyle','none');
                                temp(j:(j+k))=0;
                            end
                            hold off;
                            ylabel(['State (' units ')']);
                            xlabel(['Field Index']);
                        case 2 %2D
                            [cbar, clims] = smart_colorbar([0 1], hot(255));
                            fi_map(occupancy==0)=clims(1);
                            imagesc(centers{:},fi_map');
                            colormap(cbar);
                            set(gca,'CLim',clims);
                            
                            hold on;
                            line(...
                                [-floor(0.2*range(pos(:,1))/5)*5 0]+0.95*max(pos(:,1)),...
                                [1 1]*(min(pos(:,2)-0.1*range(pos(:,2)))),...
                                'Color','k','LineWidth',3);
                            text(0.95*max(pos(:,1)),...
                                min(pos(:,2))-0.15*range(pos(:,2)),...
                                [num2str(floor(0.2*range(pos(:,1))/5)*5) ' ' units],...
                                'VerticalAlignment','cap','HorizontalAlignment','right');
                            axis off equal;
                            set(gca,'YDir','normal');
                            hold off;
                            freezeColors;
                        case 3 %3D
                            set(gcf,'Renderer','OpenGL');
                            [cbar, clims] = smart_colorbar([0 1], hot(255));
                            fi_map(occupancy==0)=clims(1);
                            j = {0 1 2};
                            k = find(occupancy>0)';
                            view(3);
                            for k=k;
                                [j{1} j{2} j{3}] = ind2sub(size(fi_map),k);
                                clr = round((fi_map(k)-min(clims))/range(clims)*size(cbar,1))+1;
                                if clr>size(cbar,1), clr=size(cbar,1); end;
                                plotcube([1 1 1]*binside,...
                                    cellfun(@(x,y)x(y),centers,j)-[0.5 0.5 0.5]*binside,...
                                    cbar(clr,:),...
                                    min((fi_map(k)-min(clims))/range(clims),1)*0.95+0.025,...
                                    'EdgeAlpha',0,...
                                    'BackfaceCull',1);
                            end
                            l = floor(max(range(pos))*0.2/10)*10;
                            line([0 0 0;l 0 0]+min(pos(:,1))-0.01*range(pos(:,1)),...
                                [0 0 0;0 l 0]+min(pos(:,2))-0.01*range(pos(:,2)),...
                                [0 0 0;0 0 l]+min(pos(:,3))-0.01*range(pos(:,3)),...
                                'Color','k','LineWidth',3);
                            
                            text(min(pos(:,1))-0.02*range(pos(:,1)),...
                                min(pos(:,2))-0.02*range(pos(:,2)),...
                                min(pos(:,3))-0.02*range(pos(:,3)),...
                                [num2str(l) ' ' units],...
                                'VerticalAlignment','Cap','HorizontalAlignment','right');
                            hold off;
                            axis off equal;
                        otherwise % Can't plot
                            warning('pass_index:bad_dim_subplot','Cannot plot over 3d');
                    end
                end
            case 'Scatter plot'
                scatter(repmat(spk_pass_index,[2 1]),rad2deg([mod(spk_theta_phase,2*pi);mod(spk_theta_phase,2*pi)+2*pi]),'o','filled','SizeData',5);
                xlim([-1 1]);ylim([0 2*360]);
                xlabel('Pass Index');ylabel('LFP Phase (^o)');
                
                x = linspace(-1,1,750);
                phi = mod(2*pi*s*x+b,2*pi);
                k = find(abs(diff(phi))>pi);
                phi(k) = NaN;
                
                hold on;
                text(x(floor(750*0.75)),rad2deg(phi(floor(750*0.75))+2*pi)+25,['rho=' sprintf('%2.2f',rho)],'Color',[1 0 0],'FontWeight','bold','BackgroundColor','w');
                text(x(floor(750*0.75)),rad2deg(phi(floor(750*0.75))+2*pi)-25,['p=' sprintf('%2.2f',p)],'Color',[1 0 0],'FontWeight','bold','BackgroundColor','w');
                plot(x,rad2deg(phi+2*pi),'r','LineWidth',2);
                plot(x,rad2deg(phi),'r','LineWidth',2);
                
                hold off;
            case 'Density map'
                imagesc(dens_centers{1},rad2deg([dens_centers{2} dens_centers{2}+2*pi]),[density';density']);
                set(gca,'YDir','normal');
                xlabel('Pass Index');ylabel('LFP Phase (^o)');
                text(1,730,[sprintf('%3.2f',max(density(:))) 'Hz'],...
                    'HorizontalAlignment','Right',...
                    'VerticalAlignment','baseline');
                colormap jet;
                freezeColors;
            otherwise
               warning('pass_index:bad_plot','Cannot recognize plot. Skipping...');
        end
        title(plots{i});
    end
end
end

