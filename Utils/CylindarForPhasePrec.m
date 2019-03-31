function [spks_VEL4LFP,occ4Ph,fieldbound]=CylindarForPhasePrec(data_video_spk,SmoothRateMap,bound)
% CylindarForPhasePrec formats data in a open field for phase precession 
figs=0;
if figs==1
    figure;imagesc(SmoothRateMap);hold on;plot(bound(:,1),bound(:,2),'r');colormap jet
    figure;plot(data_video_spk(:,2),data_video_spk(:,3),'y');hold on
end
if length(bound)>2
    % SET UP PASSES THROUGH FIELD FOR PHASE PRECESSION
    ms=data_video_spk(:,1)*1000;
    % rescale xy to size of rate map so you can find xy cors in field
    frames=[rescale(data_video_spk(:,2),1,length(SmoothRateMap)),...
        rescale(data_video_spk(:,3),1,length(SmoothRateMap))];
    % find in durations field
    in=inpolygon(frames(:,1),frames(:,2),bound(:,1),bound(:,2));
    % pull out each pass and store in cell array
    dsig=diff([0 (abs(in')>=eps) 0]);
    startIndex=find(dsig>0);endIndex=find(dsig<0)-1;
    stringIndex=(endIndex-startIndex+1>=0);
    startIndex=startIndex(stringIndex);
    endIndex=endIndex(stringIndex);
    indices=zeros(1,max(endIndex)+1);
    indices(startIndex)=1;
    indices(endIndex+1)=indices(endIndex+1)-1;
    indices=find(cumsum(indices));
    in=zeros(length(in),1);in(indices',1)=1;
    pass=find(diff([0 indices])>1==1);
    if isempty(pass) || length(pass)==1
        spks_VEL4LFP=NaN;
        occ4Ph=NaN;
        fieldbound=NaN;
    else
        if pass(1)~=1
            pass=[1,pass];
        end
        for j=1:length(pass)
            if j+1>length(pass);break;end
            tempidx=indices(pass(j):pass(j+1)-1);
            storets{j,1}=[ms(tempidx),data_video_spk(tempidx,1),...
                data_video_spk(tempidx,2:end)];
        end
        % CRITERIA FOR PASSES THROUGH FIELD:
        %   FILTER OUT < average 3CM/SEC per pass
        %   FILTER OUT LESS THAN 200MS
        %   FILTER OUT LESS THAN 5 SPIKES OVER AT LEAST 4 THETA CYCLES
        %   MAX ISI <1000ms
        x=[];
        for j=1:length(storets)
            temp=storets{j,1};
            if figs==1
                plot(temp(:,3),temp(:,4),'k');hold on
            end
            if mean(temp(:,6))>3 && (temp(end,1)-temp(1,1))>200 ...
                    && max(diff(temp(:,1)))<1000
                if figs==1
                    plot(temp(:,3),temp(:,4),'r');hold on
                end
                % make x linear
                x=[x;linspace(0,1,size(temp,1))'];
                continue
            else
                storets{j,1}=[];
            end
        end
        
        % combine cell array of passes and create spike matrix
        % & occ matrix for use in PHprecession
        tempframes=vertcat(storets{:});
        if isempty(tempframes)
            spks_VEL4LFP=NaN;
            occ4Ph=NaN;
            fieldbound=NaN;
        else
            spks_VEL4LFP=tempframes(:,1)/1000;
            tempframes(:,3)=x;
            occ4Ph=tempframes(tempframes(:,end)==0,[1,3:end]);
            occ4Ph(:,1)=occ4Ph(:,1)/1000;
            occ4Ph(:,3)=zeros(size(occ4Ph,1),1);
            fieldbound=[0 1];
        end
    end
else
    spks_VEL4LFP=data_video_spk(data_video_spk(:,6)==1,:);
    occ4Ph=data_video_spk(data_video_spk(:,6)==0,:);
    fieldbound=NaN;
end
end

