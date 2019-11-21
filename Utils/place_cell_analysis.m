% Analysis functions for place field data.
% list of methods
% SpatialInformation
% Sparsity
% getPlaceFields
% IntraTrialStability
%
% Ryan Harvey 2019

classdef place_cell_analysis
    methods(Static)
        
        
        function spatial_information = SpatialInformation(varargin)
            % Computes the spatial information of a cell (bits/spike)
            %
            %
            % See parameters below.
            %
            % See Adrien Peyrache 2008 Methods
            %
            %  PARAMETERS
            %
            %   occupation_thresh   Bins that had less than this number of seconds of
            %                       occupation are not included in the score. (default 0)
            %
            %   n_thresh            min number of spikes before value is meaningless(default 50)
            %
            %   ratemap             pre-computed ratemaps (ncells,nbins)
            %
            %   occupancy           occupancy maps
            %
            %   n_spikes            number of spikes 

            % Modified from TSToolbox by Ryan Harvey 2019
            
            p = inputParser;
            
            p.addParameter('occupation_thresh', 0, @isnumeric);
            p.addParameter('n_thresh', 50, @isnumeric);
            p.addParameter('ratemap', [], @isnumeric)
            p.addParameter('occupancy', [], @isnumeric)
            p.addParameter('n_spikes', [], @isnumeric)
            
            p.parse(varargin{:});
            
            occupation_thresh = p.Results.occupation_thresh;
            ratemap = p.Results.ratemap;
            n_thresh = p.Results.n_thresh;
            occupancy = p.Results.occupancy;
            n_spikes = p.Results.n_spikes;
            
            % remove low occupancy bins
            ratemap(occupancy <= occupation_thresh) = NaN;
            
            % linearize
            ratemap = ratemap(:);
            occupancy = occupancy(:);
            
            % normalize to probability
            occupancy = occupancy/nansum(occupancy);
            f = nansum(occupancy.*ratemap);
            ratemap = ratemap/f;
            ix = ratemap~=0;
            SB = (occupancy(ix).*ratemap(ix)).*log2(ratemap(ix));
            spatial_information = nansum(SB);
            
            % where spiking was too low, get rid of spatial information score
            spatial_information(n_spikes<n_thresh) = NaN; 

        end
        
        
        function sparsity = Sparsity(varargin)
            % Computes the sparsity of a cell
            %
            %     Calculate single-cell sparsity index (equivalently
            %        'lifetime sparseness') as defined in Ahmed and Mehta (Trends in
            %        Neuroscience, 2009)
            %
            %   S=(1/n sum(ri))^2 / (1/n sum(ri^2))
            %
            % Where i is a spatial bin and ri is the firing rate of the cell in bin i of an environment
            % containing a total or n spatial bins. A sparsity value of 1 implies no single-cell
            % sparseness. A sparsity value approaching 0 is indicative of maximal single-cell
            % sparseness, and implies a greater amount of spatial information in each spike emitted by
            % that cell.
            %
            % See parameters below.
            %
            %
            %   PARAMETERS
            %
            %   occupation_thresh   Bins that had less than this number of seconds of
            %                       occupation are not included in the score. (0)
            %   ratemap             pre-computed ratemap
            %
            %   occupancy           occupancy map
            %
            % Ryan Harvey 2019
            
            p = inputParser;
            
            p.addParameter('occupation_thresh', 0, @isnumeric);
            p.addParameter('ratemap', [], @isnumeric)
            p.addParameter('occupancy', [], @isnumeric)
            
            p.parse(varargin{:});
            
            occupation_thresh = p.Results.occupation_thresh;
            ratemap = p.Results.ratemap;
            occupancy = p.Results.occupancy;
            
            % normalize to probability
            occupancy = occupancy / sum(occupancy(:));
            
            ratemap(occupancy <= occupation_thresh) = NaN;
            
            ratemap=ratemap(:);
            
            num = (nansum(ratemap)./ length(ratemap)).^2;
            den = nansum(ratemap.^2)./ length(ratemap);
            sparsity = num / den;
        end
        
        
        function [fields]=getPlaceFields_2d(varargin)
            % getPlaceFields_2d: locate firing fields in 2d map
            % First uses uses Kmeans and contourc (segmentImage.m - RH)
            % Then removes fields based on below params
            %
            %   Input:
            %       ratemap: symmetrical 2d map (will need updated for non-symmetrical maps)
            %       minPeakRate: min rate in hz
            %       minFieldWidth: min field width in bins
            %       maxFieldWidth: max field width in bins 
            %       maze_size_cm: size of your maze in cm
            %       debugging_fig: 0 or 1 if you want to test the code
            %
            %   Output: 
            %       fields.fieldwidth: field width in cm (largest diameter)
            %       fields.area: area of field cm^2
            %       fields.bounds: x y boundaries of each field
            %       fields.masked_field: binary mask of field location
            %       fields.peakFR: peak rate in each field
            %       fields.peakLoc: ij location of peak bin
            %       fields.com: ij center of mass (or center of field)
            %
            % Note: If no field is found, the default is to output the
            % entire rate map as a firing field. For your later place field
            % criteria, you can just specify upper and lower boundaries for
            % how large you think a place field is (>10cm & <50cm). 
            %
            % Ryan Harvey Nov. 2019
            
            p = inputParser;
            addParameter(p,'ratemap',peaks,@isnumeric)
            addParameter(p,'minPeakRate',1,@isnumeric)
            addParameter(p,'minFieldWidth',3,@isnumeric) % in bins
            addParameter(p,'maxFieldWidth',30,@isnumeric) % in bins
            addParameter(p,'maze_size_cm',100,@isnumeric)
            addParameter(p,'debugging_fig',0,@isnumeric)

            parse(p,varargin{:})
            
            ratemap = p.Results.ratemap;
            minPeakRate = p.Results.minPeakRate;
            minFieldWidth = p.Results.minFieldWidth;
            maxFieldWidth = p.Results.maxFieldWidth;
            maze_size_cm = p.Results.maze_size_cm;
            debugging_fig = p.Results.debugging_fig;

            
             if nansum(ratemap(:)) == 0
                fields.fieldwidth{1} = length(ratemap)*(maze_size_cm/length(ratemap));
                fields.area{1} = length(ratemap)*length(ratemap)*(maze_size_cm/length(ratemap));
                [r,c]=find(~isnan(ratemap));
                [k,~] = convhull(r,c);
                fields.bounds{1} = [r(k),c(k)];
                fields.masked_field{1} = ~isnan(ratemap);
                fields.peakFR{1} = max(ratemap(:));
                [r,c] = find(ratemap == fields.peakFR{1});
                fields.peakLoc{1} = [r,c];
                [x_c,y_c] = centroid(polyshape(fields.bounds{1}(:,1),fields.bounds{1}(:,2)));
                fields.com{1} = [x_c,y_c];
                fields.nfields = 1;
                return
             end
             
            if debugging_fig
                figure;
                subplot(1,3,1)
            end
            
            upscalefac = 15;
            [~,~,x,y,fieldarea,X] = segmentImage('map',ratemap,'upscalefac',upscalefac,'figs',debugging_fig);
            
            if debugging_fig
                title('kmeans & contourc')
                colorbar
            end
            
            % collect field info
            for f=1:length(x)
                %rescale coordinates
                x_temp=rescale([x{f},upscalefac+1,length(X)-upscalefac],1,(length(ratemap)));
                y_temp=rescale([y{f},upscalefac+1,length(X)-upscalefac],1,(length(ratemap)));
                x_temp = x_temp(1:end-2);
                y_temp = y_temp(1:end-2);
                
                % get largest diameter: field width cm
                fields.fieldwidth{f} = max(pdist([x_temp',y_temp']))*(maze_size_cm/length(ratemap));
                
                % get field area cm^2
                fields.area{f} = fieldarea(f);
                
                % get field boundaries
                fields.bounds{f} = [x_temp',y_temp'];
                
                % get field mask
                fields.masked_field{f} = poly2mask(x_temp,y_temp,length(ratemap),length(ratemap));
                
                % get peak rate
                fields.peakFR{f} = max(ratemap(logical(fields.masked_field{f})));
                
                % get peak location
                temp_ratemap = ratemap;
                temp_ratemap(~fields.masked_field{f}) = 0;
                if isempty(fields.peakFR{f}) || isnan(fields.peakFR{f})
                    fields.peakLoc{f} = NaN;
                    fields.peakFR{f} = NaN;
                else
                    [r,c] = find(temp_ratemap == fields.peakFR{f});
                    fields.peakLoc{f} = [r,c];
                end
                
                % get center of mass
                [x_c,y_c] = centroid(polyshape(x_temp,y_temp));
                fields.com{f} = [x_c,y_c];
            end
            
            % sort fields by firing rate
            [~,idx]=sort([fields.peakFR{:}],'descend');
            fields.fieldwidth = fields.fieldwidth(idx);
            fields.area = fields.area(idx);
            fields.bounds = fields.bounds(idx);
            fields.masked_field = fields.masked_field(idx);
            fields.peakFR = fields.peakFR(idx);
            fields.peakLoc = fields.peakLoc(idx);
            fields.com = fields.com(idx);
            
            % check fields against threshold params
            for f = 1:length(x) 
                exclude(f,1) = fields.peakFR{f} < minPeakRate | isnan(fields.peakFR{f});
            end
            for f = 1:length(x)
                exclude(f,2) = floor(fields.fieldwidth{f} / (maze_size_cm/length(ratemap))) < minFieldWidth;
            end
            for f = 1:length(x)
                exclude(f,3) = floor(fields.fieldwidth{f} / (maze_size_cm/length(ratemap))) > maxFieldWidth;
            end
            fields.fieldwidth(any(exclude,2)) = [];
            fields.area(any(exclude,2)) = [];
            fields.bounds(any(exclude,2)) = [];
            fields.masked_field(any(exclude,2)) = [];
            fields.peakFR(any(exclude,2)) = [];
            fields.peakLoc(any(exclude,2)) = [];
            fields.com(any(exclude,2)) = [];
            
            fields.nfields = length(fields.bounds);

            % if no field exist
            if isempty(fields.fieldwidth)
                fields.fieldwidth{1} = length(ratemap)*(maze_size_cm/length(ratemap));
                fields.area{1} = length(ratemap)*length(ratemap)*(maze_size_cm/length(ratemap));
                [r,c]=find(~isnan(ratemap));
                [k,~] = convhull(r,c);
                fields.bounds{1} = [r(k),c(k)];
                fields.masked_field{1} = ~isnan(ratemap);
                fields.peakFR{1} = max(ratemap(:));
                [r,c] = find(ratemap == fields.peakFR{1});
                fields.peakLoc{1} = [r,c];
                [x_c,y_c] = centroid(polyshape(fields.bounds{1}(:,1),fields.bounds{1}(:,2)));
                fields.com{1} = [x_c,y_c];
                fields.nfields = 1;
                
                if debugging_fig
                    subplot(1,3,2)
                    imAlpha=ones(size(ratemap));
                    imAlpha(isnan(ratemap))=0;
                    imagesc(ratemap,'AlphaData',imAlpha);
                    axis off
                    axis image
                    colormap(viridis(255))
                    colorbar
                    hold on
                    for f = 1:length(fields.bounds)
                        plot(fields.bounds{f}(:,1),fields.bounds{f}(:,2),'LineWidth',2)
                    end
                    title('check fields against threshold params')
                end
                return
            end

            % check for fields greater than max allowed
            for f = find([fields.fieldwidth{:}] > length(ratemap)*(maze_size_cm/length(ratemap)))
                fields.fieldwidth{f} = length(ratemap)*(maze_size_cm/length(ratemap));
                fields.area{f} = length(ratemap)*length(ratemap)*(maze_size_cm/length(ratemap));
                [r,c]=find(~isnan(ratemap));
                [k,~] = convhull(r,c);
                fields.bounds{f} = [r(k),c(k)];
                fields.masked_field{f} = ~isnan(ratemap);
                fields.peakFR{f} = max(ratemap(:));
                [r,c] = find(ratemap == fields.peakFR{f});
                fields.peakLoc{f} = [r,c];
                [x_c,y_c] = centroid(polyshape(fields.bounds{f}(:,1),fields.bounds{f}(:,2)));
                fields.com{f} = [x_c,y_c];
            end
            
            if debugging_fig
                subplot(1,3,2)
                imAlpha=ones(size(ratemap));
                imAlpha(isnan(ratemap))=0;
                imagesc(ratemap,'AlphaData',imAlpha);
                axis off
                axis image
                colormap(viridis(255))
                colorbar
                hold on
                for f = 1:length(fields.bounds)
                    plot(fields.bounds{f}(:,1),fields.bounds{f}(:,2),'LineWidth',2)
                end
                title('check fields against threshold params')
            end
            
            % remove fields with the same field boundaries and keep the one with
            % the highest peak rate and the smallest field width
            for f = 1:length(fields.peakLoc)
                fielddets(f,:) = [fields.peakLoc{f},fields.fieldwidth{f}];
            end
            [~,idx]=sort(fielddets(:,3),'ascend');

            [C,ia,~]=unique(fielddets(:,1:2),'rows','stable');
            Z=zeros(size(fielddets,1),2);
            Z(ia,:)=C;
            Z=Z(idx,:);
            
            exclude = find(all(Z==0, 2));
            fields.fieldwidth(exclude) = [];
            fields.area((exclude)) = [];
            fields.bounds((exclude)) = [];
            fields.masked_field((exclude)) = [];
            fields.peakFR((exclude)) = [];
            fields.peakLoc((exclude)) = [];
            fields.com((exclude)) = [];
            
            if debugging_fig
                subplot(1,3,3)
                imAlpha=ones(size(ratemap));
                imAlpha(isnan(ratemap))=0;
                imagesc(ratemap,'AlphaData',imAlpha);
                axis off
                axis image
                colormap(viridis(255))
                colorbar
                hold on
                for f = 1:length(fields.bounds)
                    plot(fields.bounds{f}(:,1),fields.bounds{f}(:,2),'LineWidth',2)
                end
                title('remove fields with the same field boundaries')
            end
            
            polyvec = [];
            for f = 1:length(fields.bounds)
                polyvec = [polyvec, polyshape(fields.bounds{f}(:,1),fields.bounds{f}(:,2))];
            end

            TF = overlaps(polyvec);
            if sum(TF) == 1
                fields.nfields = length(fields.bounds);
                return
            end
            exclude = [];
            for r = 1:length(TF)
               for c = 2:length(TF) 
                    if TF(r,c) == 1
                        polyout = intersect([polyshape(fields.bounds{r}(:,1),fields.bounds{r}(:,2)),...
                            polyshape(fields.bounds{c}(:,1),fields.bounds{c}(:,2))]);
                        if round(polyarea(polyout.Vertices(:,1),polyout.Vertices(:,2)),2)==...
                                round(polyarea(fields.bounds{c}(:,1),fields.bounds{c}(:,2)),2)
                            exclude(r,c) = 1;
                        end
                    end
               end
            end
            for f = length(exclude):-1:1
                if exclude(1,f) == 1
                    exclude_(f) = 1;
                end
            end
            
            fields.fieldwidth(logical(exclude_)) = [];
            fields.area(logical(exclude_)) = [];
            fields.bounds(logical(exclude_)) = [];
            fields.masked_field(logical(exclude_)) = [];
            fields.peakFR(logical(exclude_)) = [];
            fields.peakLoc(logical(exclude_)) = [];
            fields.com(logical(exclude_)) = [];
            
            fields.nfields = length(fields.bounds);

            if debugging_fig
                subplot(1,3,3)
                imAlpha=ones(size(ratemap));
                imAlpha(isnan(ratemap))=0;
                imagesc(ratemap,'AlphaData',imAlpha);
                axis off
                axis image
                colormap(viridis(255))
                colorbar
                hold on
                for f = 1:length(fields.bounds)
                    plot(fields.bounds{f}(:,1),fields.bounds{f}(:,2),'LineWidth',2)
                end
                title('remove fields with the same field boundaries')
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
            addParameter(p,'ratemap',@isnumeric)
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
            

            stdRates = std(ratemap,[],2)';
            ratemap = mean(ratemap,2)';

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
                % remove field peaks with a half standard dev higher than the mean
                % (unreliable fields)
                exclude = [exclude; find( ratemap(i,locs) <  stdRates(i,locs)*.5)'];
                
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
        
        
        function r=IntraTrialStability(mat,track,track_length)
            
            if isempty(mat)
                r=NaN;
                return
            end
            
            if track==1
                nBinsy=1;
                filtWidth=[1,5];
            elseif track==0
                nBinsy=round(track_length/3);
                filtWidth = [5 5];
            end
            
            
            split=round(size(mat,1)/2);
            S{1}=mat(1:split,:);
            S{2}=mat(split+1:end,:);
            MinY=min(mat(:,3));MaxY=max(mat(:,3));MinX=min(mat(:,2));MaxX=max(mat(:,2));
            nBinsx = round(track_length/2);
            edges{1}=linspace(MinY,MaxY,nBinsy+1);
            edges{2}=linspace(MinX,MaxX,nBinsx+1);
            
            for i=1:2
                mat=S{i};
                Omatrix = hist3([mat(mat(:,6)==0,3),mat(mat(:,6)==0,2)],'Edges',edges);
                Omatrix(end,:) = [];
                Omatrix(:,end) = [];
                occ = Omatrix/30;
                Smatrix = hist3([mat(mat(:,6)==1,3),mat(mat(:,6)==1,2)],'Edges',edges);
                Smatrix(end,:)=[];
                Smatrix(:,end)=[];
                FilledRateMatrix=Smatrix./occ;
                FilledRateMatrix(isinf(FilledRateMatrix))=0;
                imageFilter=fspecial('gaussian',filtWidth,1);
                map{i}=nanconv(FilledRateMatrix,imageFilter, 'nanout');
            end
            map1=map{1};
            map2=map{2};
            map1(isnan(map1))=0;
            map2(isnan(map2))=0;
            r=corrcoef(map1,map2);
            r=r(2,1);
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
                        kempter_lincirc(data.position.x(data.position.lap==i),...
                        rad2deg(data.position.phase(data.position.lap==i)));
                else
                    circ_lin_corr(i,1)=NaN;pval(i,1)=NaN;slope_deg(i,1)=NaN;phi0_deg(i,1)=NaN;RR(i,1)=NaN;
                end
            end
            ThPrecess.lapCorrelation=nanmean(circ_lin_corr);
            ThPrecess.lapPhaseOffset=nanmean(phi0_deg);
            
            [ThPrecess.circLinCorr,ThPrecess.pval,ThPrecess.slopeCpU,ThPrecess.phaseOffset]=kempter_lincirc(data.position.x,data.position.phase);
        end
    end
end