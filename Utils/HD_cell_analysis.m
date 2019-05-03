classdef HD_cell_analysis
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods(Static)
        
        %Distributive Ratio
        function DR = distributiveRatio(ratemap,frames,hdTuning,mazesize)
            
            % thetabins=linspace(0,360,60);
            thetabins=0:6:360;
            
            nBinsx = round(mazesize/3); nBinsy = round(mazesize/3);
            MinY = min(frames(:,3));
            MaxY = max(frames(:,3));
            MinX = min(frames(:,2));
            MaxX = max(frames(:,2));
            edges{1} = linspace(MinY, MaxY, nBinsy+1);
            edges{2} = linspace(MinX, MaxX, nBinsx+1);
            
            % Build
            for i=1:length(thetabins)-1
                % find frames in specific theta bin
                frames_by_theta=frames(frames(:,4)>thetabins(i) & frames(:,4)>thetabins(i+1),:);
                % bin xy
                theta_map = hist3([frames_by_theta(:,3) frames_by_theta(:,2)],'Edges',edges);
                theta_map(end,:) = [];
                theta_map(:,end) = [];
                
                PR(i)=nansum(nansum(ratemap.*theta_map))/nansum(theta_map(:)); %Calculates predicted rate
            end
            
            for i=1:length(hdTuning)
                tempDR(i)=log((1+hdTuning(i)))/(1+PR(i)); 
            end
            
            DR=nansum(abs(tempDR))/length(hdTuning);
            
        end
        
        %Directional Information Content 
        %Inputs: 
        %Output: Directional information content. 
        function DIC  = computeDIC(histAng,hdTuning,OverallFR)
            %Computes the directional information content on unsmoothed ratemap for
            %directional data. Adapted from Langston et al 2010. by LB March 2018
            %
            probOrient = histAng./sum(histAng);
            reIC = hdTuning./OverallFR;
            log_IC = log2(reIC); log_IC(isinf(log_IC)) = 0;
            ICi = probOrient.*reIC.*log_IC;
            DIC = sum(ICi);
        end
        
        function stability=hd_stability(data_video_spk)
            da=pi/30;
            angBins=rad2deg(da/2:da:2*pi-da/2);
            bin_centers=movmedian(angBins,2);
            bin_centers(1)=[];
            i=1;
            ii=2;
            chunk=[];
            hdTuning=[];
            heading=data_video_spk(:,4);
%             figure;
%             plot(data_video_spk(:,2),data_video_spk(:,3),'.k')
            
            while ii<length(heading)
                histAng=histcounts(heading(i:ii),angBins);
                if sum(histAng>0)==length(angBins)-1
                    temp=data_video_spk(i:ii,:);
                    spk_a=temp(temp(:,6)==1,4);
                    a=temp(temp(:,6)==0,4);
                    [~,~,~,~,~,hdtuning]=tuningcurve(a,spk_a,30);
                    hdTuning=[hdTuning;hdtuning];
                    chunk=[chunk;i,ii];
%                     hold on
%                     plot(data_video_spk(i:ii,2),data_video_spk(i:ii,3))
                    i=ii+1;
                    ii=ii+1;
                end
                ii=ii+1;
            end
           
            store=[];
            corrMat=corr(hdTuning');
             for k=1:size(corrMat,2)-1
                 data=diag(corrMat,k);
                 store=[store;data];
             end
            
             stability=nanmean(store);
             
        end
        
    end
end

