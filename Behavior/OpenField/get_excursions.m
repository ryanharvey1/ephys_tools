
cd('D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data')
load('params_V17.mat') %%loads table created from OF_preprocess.
clearvars -except params
figures=1;
fr=30;
param_idx=params.subID; %serves as index for getParam
%Loop through subjects
for j=1:size(param_idx,1)
    
    %Grab variables for easier use
    back=params.backCM{j}; %using center of back xy coordinates
    ts=params.ts{j};
    
    excursion{1}=[];
    
    stopIdx=contiguousframes(params.pathIV{j}<3,30); %find segments of coordinates where rat is moving less than 3cm/s
    [startStop,endStop,~]=findgroups(params.stopIdx{j}); %find the ends of the stops to determine begining of movements
    idxEnd=zeros(size(back,1),1); idxStart=zeros(size(back,1),1); %initialize the logical index. 
    idxEnd(endStop',1)=1; %populate the logical index to find the row corresponding to end of stops. 
    idxStart(startStop',1)=1; %populate the logical index to find the row corresponding to end of stops. 
    
    stopIdx(end+1)=stopIdx(end);
    % loop through home bases
    for hb=1:size(params.HBBound{j},2)
        stopIn = inpolygon(back(:,1),back(:,2),...
            params.HBBound{j}{1, hb}(:,1),params.HBBound{j}{1, hb}(:,2)); %find which end of stops were in the home base

        % cross reference slow-down times with times in the home base
        idx=stopIn & stopIdx;
        [startStop,endStop,~]=findgroups(~idx);
        
        tempHBX=sin(0:2*pi/1000:2*pi)*15+params.HBcenter{j}{1,hb}(:,1); 
        tempHBY=cos(0:2*pi/1000:2*pi)*15+params.HBcenter{j}{1,hb}(:,2); 
        
        % loop through each segment and keep if: 
        %   1. most of the time is not spent in the homebase
        %   2. the start of the segment is in the home base
        %   3. the end is in the home base
        for i=1:length(startStop)
            x=back(startStop(i):endStop(i),1);
            y=back(startStop(i):endStop(i),2);
            
            if size(x,1) < 2 || size(y,1) < 2 %skip segments that are only one point
                continue 
            end
            
            [ ~, pathIV, ~] = compute_OFpathCalc(x,y,fr);
            
            if max(pathIV) < 3 %skip segments where the rat doesn't move
                continue
            end
            
            if sum(inpolygon(x,y,params.HBBound{j}{1, hb}(:,1),params.HBBound{j}{1, hb}(:,2)))<length(x)*.9 &&...
                    inpolygon(x(1),y(1),tempHBX,tempHBY) &&... %params.HBBound{j}{1, hb}(:,1),params.HBBound{j}{1, hb}(:,2)
                    inpolygon(x(end),y(end),tempHBX,tempHBY)
                keep_seg{i}=[x,y];
            end

        end
        % remove empty cells
        keep_seg=keep_seg(~cellfun('isempty',keep_seg));
        
        % plot for sanity check
        if figures==1
            figure
            plot(back(:,1),back(:,2),'color',[.7 .7 .7]);hold on
            xlabel('x')
            ylabel('y')
            plot(params.HBBound{j}{1, hb}(:,1),params.HBBound{j}{1, hb}(:,2),'g')
            
            for i=1:length(keep_seg)
                pause(2)
                plot(keep_seg{i}(:,1),keep_seg{i}(:,2),'linewidth',2)
                scatter(keep_seg{i}(1,1),keep_seg{i}(1,2),'g')
                scatter(keep_seg{i}(end,1),keep_seg{i}(end,2),'r')
            end
            darkBackground(gcf,[0.2 0.2 0.2],[1 1 1])
        end
        
        seg{hb}=keep_seg;
        clear keep_seg
        
    end
end