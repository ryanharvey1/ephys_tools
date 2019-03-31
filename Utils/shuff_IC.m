% shuff_IC

% load('data.mat')
close all;clc
rats={'RH13','RH14','LS21','LS23','LE2821','LE2823','RH11','RH16','LS17','LS19','LE2813'};
% figure;box off
IC_SHUFF=[];
for r=1:length(rats)
    sessions=fieldnames(data.(rats{r}));
    for s=1:length(sessions)
        % extract events and frames
        events=data.(rats{r}).(sessions{s}).events;
        frames=data.(rats{r}).(sessions{s}).frames;
        
        % restrict to linear track session
        if sum(events==[1;1])~=2 % if more than one session
            frames=frames(frames(:,1)>events(1,1) & frames(:,1)<events(2,1),:);
        end
        
        % find track length
        date=char(sessions{s});
        date=[date(2:5),'-',date(6:7),'-',date(8:9),'_',date(10:11),'-',date(12:13),'-',date(14:15)];
        track_length=TrackLength(['/Volumes/Ryan_4TB/Place_Cell_Data/RawPAE_PlaceCell',filesep,rats{r},filesep,date]);
        
        disp([rats{r},filesep,date])
        
        % calc head angle
        ExtractedAngle=XYangle(frames(:,2),frames(:,3));
        
        % calc vel
        vel_cmPerSec=abs(diff(frames(:,2)))*track_length/(range(frames(:,2)))*30;
        
        frames(:,[4:5])=[[ExtractedAngle;ExtractedAngle(end)],[vel_cmPerSec;vel_cmPerSec(end)]];
        
        % extract spike times
        S=data.(rats{r}).(sessions{s}).Spikes;
        
        for i=1:length(S)
            
            SpikeFile=S{i};
            
            if sum(events==[1;1])~=2 % if more than one session
                SpikeFile=SpikeFile(SpikeFile(:,1)>events(1,1) & SpikeFile(:,1)<events(2,1),:);
            end
            
            if length(SpikeFile)<50
                continue
            end
            
            in=contiguousframes(frames(:,5)<2,60);
            data_video_nospk=[frames,in];
            
            % INTERPOLATE SPIKES TO TIMESTAMPS, POSITION, AND VEL DATA
            X=interp1(data_video_nospk(:,1),data_video_nospk(:,2),SpikeFile,'linear');
            Y=interp1(data_video_nospk(:,1),data_video_nospk(:,3),SpikeFile,'linear');
            A=interp1(data_video_nospk(:,1),data_video_nospk(:,4),SpikeFile,'linear');
            VEL=interp1(data_video_nospk(:,1),data_video_nospk(:,5),SpikeFile,'linear');
            VELidx=interp1(data_video_nospk(:,1),data_video_nospk(:,6),SpikeFile,'linear');
            
            % CONCAT AND SORT
            data_video_spk=sortrows([[SpikeFile X Y A VEL VELidx ones(size(SpikeFile,1),1)];[data_video_nospk,zeros(length(data_video_nospk),1)]],1);
            
            % VELO FILTER BASED ON INDEX CREATED ABOVE
            data_video_nospk(logical(in),:)=[];
            data_video_nospk(:,6)=[];
            data_video_spk(data_video_spk(:,6)==1,:)=[];
            data_video_spk(:,6)=[];

            [right,left,DirectionalityIndex,Displacement,nlaps]=RightVsLeft(data_video_nospk,data_video_spk,track_length,30);
            secIdx=data_video_nospk(30:30:length(data_video_nospk),1);
            
            IC=SHUFF(right,secIdx,track_length);
            IC_SHUFF=[IC_SHUFF;IC];
            IC=SHUFF(left,secIdx,track_length);
            IC_SHUFF=[IC_SHUFF;IC];          
        end
    end
end
NintyFith_Percentile = prctile(IC_SHUFF,95)
% plot
h = histogram(IC_SHUFF,'FaceColor','k');hold on;
plot([prctile(IC_SHUFF,95);prctile(IC_SHUFF,95)],[0;max(h.BinCounts)])
hold off;

function IC=SHUFF(structure,secIdx,track_length)
data_video_spk=structure.dataspks;
[I,row]=ismember(data_video_spk(:,1),secIdx);
I=find(I);
I=[1;I;length(data_video_spk)];

for s=1:length(I)
    if s==length(I)
        datacell{s,1}=data_video_spk(I(s):I(end),:);
        break
    end
    datacell{s,1}=data_video_spk(I(s):I(s+1)-1,:);
end

tempframes=data_video_spk;
shuff_max=(length(secIdx)-20);
bispk=[];
IC=[];
for ishuff=1:400
    shift=randi([20 round(shuff_max)]);
    tempcell=circshift(datacell,shift);
    for open=1:length(tempcell)
        bispk=[bispk;tempcell{open}(:,6)];
    end
    tempframes(:,6)=bispk;
    bispk=[];
    
    [ SmoothRateMap,nBinsx,nBinsy,occ,~] = bindata(tempframes(tempframes(:,6)==0,:),30,tempframes(tempframes(:,6) == 1,:),'yes',track_length);
    rY=reshape(SmoothRateMap,nBinsx*nBinsy,1);
    rY(isnan(rY) & isinf(rY))=0;
    occRSHP=reshape(occ,nBinsx*nBinsy,1);
    pX=occRSHP./sum(occRSHP);
    [nBins,nCells]=size(rY);
    relR=rY./kron(ones(nBins,1),pX'*rY);
    log_relR=log2(relR);
    log_relR(isinf(log_relR))=0;
    IC = [IC; sum(kron(pX,ones(1,nCells)).*relR.*log_relR)];
end
end

% function IC=SHUFF(vid,spk,track_length)
% spksf=spk(:,6);
% IC=[];
% binary=zeros(length(spksf),400);
% for i=1:400
%     binary(:,i)=circshift(spksf,randi(length(spk),1));
% end
% tic
% parfor x = 1:400
% %     spk(:,6)=circshift(spksf,randi(length(spk),1));
%     [ SmoothRateMap,nBinsx,nBinsy,occ,~] = bindata(vid,30,spk(logical(binary(:,x)),:),'yes',track_length);
%     rY=reshape(SmoothRateMap,nBinsx*nBinsy,1);
%     rY(isnan(rY))=0;rY(isinf(rY))=0;
%     occRSHP=reshape(occ,nBinsx*nBinsy,1);occSUM=sum(occRSHP);pX=occRSHP./occSUM;
%     [nBins,nCells]=size(rY);relR=rY./kron(ones(nBins,1),pX'*rY);
%     log_relR=log2(relR);log_relR(isinf(log_relR))=0;
%     InformationContent=sum(kron(pX,ones(1,nCells)).*relR.*log_relR);
%
% %     IC = [IC; InformationContent];
% end
% toc
% end
