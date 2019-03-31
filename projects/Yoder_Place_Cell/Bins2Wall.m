%% Bins2Wall
% THIS SCRIPT READS IN RYAN YODER'S RATEMAP MATRIX DATA AND OUTPUTS BINS TO WALL AND OTHER STATS
%
% Input:
%       - path to data output from Ryan Yoder's programs
%
% Output:
%       - rate map
%       - rate map stats
%
% Created by Ben C June 2015
% Updated by Ryan H Nov. 2016

clear, close, clc

% ##############
ratemap=0;
matfile=0;
path='/Users/ryanharvey/Downloads/Place Cells - Tilted Mice'; % PATH TO PARENT FOLDER
% ##############

% ADD PATHS TO TOOLBOXES
addpath(genpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/NSMA_Toolbox')); % add NSMA to path for FindFiles
addpath(genpath('/Users/RyanHarvey/Dropbox/MATLAB/NSMA_Toolbox'));
addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\Cell analysis'));
addpath(genpath('/Users/ryanharvey/Dropbox/MATLAB/Cell analysis'));
addpath(genpath('/Users/ryanharvey/Dropbox/MATLAB/chronux_2_11'));
addpath(genpath('/Users/ryanharvey/Dropbox/MATLAB/spikeCode'));
addpath(genpath('/Users/ryanharvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis'));

parent=dir(path); row=1; counter = 1;

for k=1:length(parent)
    if parent(k).isdir && parent(k).name(1) ~= '.' % IF DIR NOT '.'
        cd([path filesep parent(k).name filesep 'ASCII']);
        
%         txtFiles=FindFiles('*txtrmap*.txt');
        txtFiles=dir('*txtrmap*.txt');
        txtFiles={txtFiles.name}';

        
        for kk=1:length(txtFiles)
            
            % load rate map array
            data = importdata(txtFiles{kk});
            
            % extract coords, spikes, angle, and direction from ReadData output
            RateMap = data.data; % rate map
            
            % LOCATE BORDERS
            [rows,columns]=size(RateMap);
            for i=1:rows
                for J=1:columns
                    if isnan(RateMap(i,J))==0
                        index1(i)=i;
                    end
                end
            end
            index1=index1(index1~=0); y_max=index1(1); y_min = numel(RateMap(:,1));
            
            for i=1:columns
                for J=1:rows
                    if isnan(RateMap(J,i))==0
                        index(i)=i;
                    end
                end
            end
            index=index(index~=0); x_max=index(1); x_min = numel(RateMap(1,:));
            
            
            %             [ sparsity,Coherence,InformationContent,Peak_Rate] = altBins2Wall(RateMap,x_min,y_min);
            
            
            % calculate rate map statistics
            r = reshape(RateMap,x_min*y_min,1);
            r(isnan(r))=[];
            Peak_Rate = max(r);
            
            % Removes peak rates that are too high for biology
            if numel(find(r > (Peak_Rate*0.20)))<=1 % find all spikes >20% of max firing rate
                r(r==Peak_Rate)=[]; RateMap(RateMap==Peak_Rate)=NaN; Peak_Rate = max(r);
            end
            
            % SMOOTH RATE MAP
            Ratemap1=RateMap(y_max:y_min,x_max:x_min);
            Ratemap1(isnan(Ratemap1))=0;
            
            %             [SmoothRateMap] = smooth(Ratemap1,1,5,1,5);
            SmoothRateMap=Ratemap1;
            
            % RESHAPE
            
            [rowz,col]=size(SmoothRateMap);
            r2 = reshape(SmoothRateMap,rowz*col,1);
            
            % SPARSITY FROM RESHAPED SMOOTHED RATE MAP
            [~,~,sparsity]=Sparsity(r2);
            
            % INFORMATION CONTENT
            [InformationContent, infomat] = InformationPerSpike({r2},r2./sum(r2)); % from NSMA toolbox
            
            % COHERENCE
            if length(r2)>1
                Coherence=corr2(Ratemap1,SmoothRateMap);
            else
                Coherence=NaN;
            end
            
            %             figure(1);pcolor(Ratemap1);figure(2);pcolor(SmoothRateMap);figure(3);pcolor(RateMap);
            
            % AVERAGE RATE
            AVGRate=mean(r2);
            
            if isempty(Peak_Rate)==1
                Peak_Rate=NaN; % turning [] to NaN for table
            end
            % 20 percent % 20 % of max firing rate % find all spikes >20% of max firing rate
            r_20perc = Peak_Rate*0.20; r_20percAll = find(r > r_20perc); NumbActiveBins20 = numel(r_20percAll);
            
            if NumbActiveBins20<1
                NumbActiveBins20=NaN;
            end
            AllBins = numel(r); % number of bins % Field Width % percent of active bins
            FieldWidth20 = NumbActiveBins20; PercActiveBins20 = (NumbActiveBins20/AllBins)*100;
            
            % 50 percent % 50% of max firing rate % find all spikes >50% of max firing rate
            r_50perc = Peak_Rate*0.50; r_50percAll = find(r > r_50perc); NumbActiveBins50 = numel(r_50percAll);
            
            if NumbActiveBins50<1
                NumbActiveBins50=NaN;
            end
            FieldWidth50 = NumbActiveBins50; PercActiveBins50 = (NumbActiveBins50/AllBins)*100;
            
            if NumbActiveBins20>=1
                [y_max_rate, x_max_rate] = find(RateMap == Peak_Rate); % xy position of peak firing spike
                x_top = x_min - x_max_rate; x_bottom = x_max_rate - x_min;
                y_top = y_min - y_max_rate; y_bottom = y_max_rate - y_min;
                distanceAll = [x_top x_bottom y_top y_bottom];
                distanceMin = min(distanceAll); distanceBins = distanceMin;
            end
            % plot raw rate map
            if ratemap==1
                fig1=figure(kk); h = pcolor(RateMap);
                colormap(jet);
                axis square tight
                hold on
                set(h, 'EdgeColor', 'none');
                colorbar;
                box off
                if NumbActiveBins20>=1
                    scatter(x_max_rate+0.5, y_max_rate+0.5, 100, 'k', 's', 'LineWidth', 1.5);
                end
                set(gca,'YDir','reverse');
                xlim([x_max-2 x_min+2]); ylim([y_max-2 y_min+2]);
            end
            
            % Plot Circle
            x=median([x_max,x_min]); y=median([y_max,y_min]); rad=(max(x_min-x_max,y_min-y_max))/2;
            th = 0:pi/179.5:2*pi; % 0 to 2*pi(6.28318530717959) at 0.0175 increments to equal 360 points
            xunit = rad * cos(th) + x; yunit = rad * sin(th) + y;
            if ratemap==1
                figure(kk)
                hold on
                h = plot(xunit, yunit);
                hold off
            end
            
            % Calc bins from peak to wall
            peak_bins2wall=NaN;
            if NumbActiveBins20>=1
                
                for ii=1:length(xunit)
                    d(ii) = sqrt((x_max_rate(1)-xunit(ii))^2+(y_max_rate(1)-yunit(ii))^2);
                end
                peak_bins2wall=min(d);
            end
            % CALCULATE BINS FROM FIELD TO WALL (20% & 50% THRESHOLDS)
            field_bins2wall20=NaN; field_bins2wall50=NaN;
            
            if NumbActiveBins20>=1
                % 20 percent
                [y_20,x_20] = find(RateMap > r_20perc);
                for j=1:length(x_20) % calculates min dist from every active bin to the walls
                    for jj=1:length(xunit)
                        d2(jj) = sqrt((x_20(j)-xunit(jj))^2+(y_20(j)-yunit(jj))^2);
                    end
                    d3(j)=min(d2);
                end
                field_bins2wall20=mean(d3);
                
                clear y_20 x_20 j jj d2 d3
                
                % 50 percent
                [y_20,x_20] = find(RateMap > r_50perc);
                for j=1:length(x_20) % calculates min dist from every active bin to the walls
                    for jj=1:length(xunit)
                        d2(jj) = sqrt((x_20(j)-xunit(jj))^2+(y_20(j)-yunit(jj))^2);
                    end
                    d3(j)=min(d2);
                end
                field_bins2wall50=mean(d3);
                
                % PLOT MIN LINE OF PEAK RATE
                if ratemap==1
                    for iii=1:length(xunit)
                        dis(iii) = sqrt((x_max_rate(1)-xunit(iii))^2+(y_max_rate(1)-yunit(iii))^2);
                    end
                    [~,I]=min(dis); % use min indexing
                    figure(kk)
                    hold on
                    plot([x_max_rate(1) xunit(I)], [y_max_rate(1) yunit(I)],'Color','k','LineWidth',6); % to plot min line using indexed min
                end
            end
            % BORDER SCORE
            Map4border=RateMap(y_max:y_min,x_max:x_min);
            [field,fieldwidth]=FindFF2D(Map4border);
            if isnan(field)==0
                Map4bordertemp=Map4border;
                Map4bordertemp(~field)=NaN;
                infieldFR=nanmean(reshape(Map4bordertemp,[],1));
%                 IF=nansum(reshape(Map4bordertemp,[],1));
                
                Map4bordertemp=Map4border;
                Map4bordertemp(logical(field))=NaN;
                outfieldFR=nanmean(reshape(Map4bordertemp,[],1));
%                 OF=nansum(reshape(Map4bordertemp,[],1));
                
%                 SNR=(IF-OF)/(IF+OF);
                
                nanPlacement=isnan(Map4border);
                Map4border(~field)=0;
                Map4border(nanPlacement)=NaN;
            else
                infieldFR=NaN;
                outfieldFR=NaN;
%                 SNR=NaN;
            end
            border_score=BorderScore(Map4border);
            
            
            [filepath, filename] = fileparts(txtFiles{kk});
            
            % FIELD DISPLACEMENT FROM SESSION 1 TO 2
            sess=char(filename);sess=sess(end);
            if sess=='1' || sess=='2'
                if sess=='1'
                    SaveField=field;
                elseif sess=='2'
                    [D,M]=Displacement(SaveField(:,:,1),field);
                else
                    clear SaveField
                end
            end
            
            % SAVE RATEMAP
            if ratemap==1
                saveas(fig1,[filepath filesep filename '_ratemap.jpg']);
            end
            
            % COMPILE STATS
            if ~exist('D','var');D=NaN;M=NaN;end
            
            Bins2wall_Stats(row,:)=[strrep(txtFiles{kk}, ',', ''),num2cell(field_bins2wall20),num2cell(field_bins2wall50)...
                num2cell(peak_bins2wall),num2cell(Peak_Rate),num2cell(AVGRate),num2cell(PercActiveBins20),num2cell(PercActiveBins50),...
                num2cell(FieldWidth20),num2cell(FieldWidth50),num2cell(sparsity),num2cell(InformationContent),num2cell(Coherence),num2cell(border_score),...
                fieldwidth,D,M,infieldFR,outfieldFR];
            row=row+1;
            
            % SAVE MAT FILE
            if matfile==1
                save([filepath filesep filename '_Data.mat']);
            end
            clearvars ('-except', 'path', 'parent', 'k', 'kk', 'txtFiles', 'row','Bins2wall_Stats','ratemap','matfile','counter','SaveField');
            close all
        end
        % SAVE STATS TABLE TO PARENT FOLDER
        [filepath, filename] = fileparts(path);
        
        Bins2wall_Stats_table=cell2table(Bins2wall_Stats,'VariableNames',...
            {'Filename','Field_to_Wall_20per','Field_to_Wall_50per','Peak_Bin_2_Wall','Peak_Rate','AVG_Rate',...
            'Percent_Active_Bins_20','Percent_Active_Bins_50','FieldWidth20','FieldWidth50','Sparsity','Info_Content',...
            'Coherence','BorderScore','fieldwidth','Displacement','DisplacementCorrelation','infieldFR','outfieldFR'});
        if ismac==1
            writetable(Bins2wall_Stats_table,[filepath filesep 'Field_Stats' filename '.csv'],'WriteRowNames',true);
        else
            writetable(Bins2wall_Stats_table,[filepath filesep 'Field_Stats' filename '.xlsx']);
        end
        disp(['Just finished iteration #', num2str(counter),' of ',num2str(length(parent))]);
        counter = counter + 1;
    end
end


