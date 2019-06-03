function shuff_pass=shuff(groupid,varargin)
% shuff: general shuffler
%
% Input:
%       groupid: cell array with same format as below (postprocessed .mat file,
%                   tetrode id, & cell number)
%
%     {'RH13_S20160808101427.mat','TT2.mat','2';
%     'RH13_S20160808101427.mat','TT2.mat','4';
%     'RH13_S20160808101427.mat','TT4.mat','5';
%     'RH13_S20160808101427.mat','TT5.mat','2';
%     'RH13_S20160808101427.mat','TT8.mat','3';
%     'RH13_S20160808103145.mat','TT4.mat','1'}
%       ...
%
%       options: ---see defaults below
%               feature:    cell array of features you want to shuffle
%                           currently supports ic, mvl, dic
%                           to add more, just add another sub-function
%               percentile: percentile threshold
%               session:    session number
%               plotting:   want plots? 1 or 0
%               nshuffle:   number of shuffles
% output:
%       shuff_pass: binary of
%
% ryan harvey 2019

% Output:
p = inputParser;
p.addParameter('feature',{'ic'});
p.addParameter('percentile',99);
p.addParameter('session',1);
p.addParameter('plotting',0);
p.addParameter('nshuffle',200);
p.addParameter('runningdir',[]);

p.parse(varargin{:});

feature = p.Results.feature;
percentile = p.Results.percentile;
session = p.Results.session;
plotting = p.Results.plotting;
nshuffle = p.Results.nshuffle;
runningdir = p.Results.runningdir;

% sessions=unique(groupid(:,1),'stable');
sessions=groupid;
shuff_pass=[];

fields_to_load={'frames','events','Spikes','spikesID',...
    'maze_size_cm','sessionID','samplerate'};

if ~isempty(runningdir)
    fields_to_load={'frames','events','Spikes','spikesID',...
        'maze_size_cm','sessionID','samplerate','linear_track'};
    dirs={'right','left'};
    
    %  [U_CA,I,J]=uniqueRowsCA([groupid(:,1),num2cell(runningdir)])
    %  U_CA(J)
    % groupid(J,1)
    % unique([groupid(:,1),num2cell(runningdir)])
    
    % sessions = cellfun(@(S) S(1:end-1),unique(strcat(groupid(:,1),...
    %     num2str(runningdir)),'stable'), 'Uniform', 0);
end
for i=1:size(sessions,1)
    
    data=load(sessions{i},fields_to_load{:});
    
    %     if isempty(runningdir)
    
    %         idx=contains(groupid(:,1),sessions{i});
    %
    %         cells_to_find=strcat(groupid(idx,2),num2str(str2double(groupid(idx,3))));
    %
    %         cell_list=strcat(data.spikesID.TetrodeNum,num2str(data.spikesID.CellNum));
    %
    %         cells=find(ismember(cell_list,cells_to_find))';
    cells_to_find=strcat(groupid(i,2),num2str(str2double(groupid(i,3))));
    cell_list=strcat(data.spikesID.TetrodeNum,num2str(data.spikesID.CellNum));
    cells=find(ismember(strrep(cell_list,' ',''),cells_to_find))';
    %     elseif ~isempty(runningdir)
    %         idx1=contains(groupid(:,1),sessions{i}) & runningdir==1;
    %         idx2=contains(groupid(:,1),sessions{i}) & runningdir==2;
    %
    %         cells=[find(ismember(strcat(data.spikesID.TetrodeNum,num2str(data.spikesID.CellNum)),...
    %             strcat(groupid(idx1,2),num2str(str2double(groupid(idx1,3))))))',...
    %             find(ismember(strcat(data.spikesID.TetrodeNum,num2str(data.spikesID.CellNum)),...
    %             strcat(groupid(idx2,2),num2str(str2double(groupid(idx2,3))))))'];
    %
    %        temprunningdir=runningdir(contains(groupid(:,1),sessions{i}));
    %     end
    %     for icell=cells
    if isempty(runningdir)
        [data_video_spk,data_video_nospk]=createframes_w_spikebinary(data,session,cells);

    elseif ~isempty(runningdir)
        data_video_spk=data.linear_track{1,1}.(dirs{runningdir(i)}){1,cells}.dataspks;
        data_video_nospk=data_video_spk(data_video_spk(:,6)==0,:);
    end
    
    for f=1:length(feature)
        truth{f,1}=feval(feature{f},data_video_spk,data,session);
    end
    
    secIdx=data_video_nospk(data.samplerate:data.samplerate:length(data_video_nospk),1);
    
    [I,~]=ismember(data_video_spk(:,1),secIdx);
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
    for ishuff=1:nshuffle
        shift=randi([20 round(shuff_max)]);
        tempcell=circshift(datacell,shift);
        for open=1:length(tempcell)
            bispk=[bispk;tempcell{open}(:,6)];
        end
        tempframes(:,6)=bispk;
        bispk=[];
        
        for f=1:length(feature)
            shuffled{f,ishuff}=feval(feature{f},tempframes,data,session);
        end
    end
    
    clear datacell
    
    shuff_pass(i,:)=cell2mat(truth)'>prctile(cell2mat(shuffled)',percentile);
    
    disp([num2str((length(shuff_pass)/length(groupid(:,1)))*100) ' percent done'])
end
end
% end

function IC=ic(tempframes,data,session)
if isfield(data,'linear_track')
    [ratemap,~,~,occ,~]=bindata(tempframes(tempframes(:,6)==0,:),...
        data.samplerate,tempframes(tempframes(:,6)==1,:),'yes',data.maze_size_cm(session));
else
    [ratemap,~,~,occ,~]=bindata(tempframes(tempframes(:,6)==0,:),...
        data.samplerate,tempframes(tempframes(:,6)==1,:),'no',data.maze_size_cm(session));
end
IC=place_cell_analysis.SpatialInformation('ratemap',...
    ratemap,'occupancy',occ,'n_spikes',sum(tempframes(:,6)));
end
function r=mvl(tempframes,data,session)

[r,I,Ispk,peakrate,prefdirec,hdTuning]=tuningcurve(tempframes(tempframes(:,6)==0,4),...
    tempframes(tempframes(:,6)==1,4),data.samplerate);
end
function Ispk=dic(tempframes,data,session)

[r,I,Ispk,peakrate,prefdirec,hdTuning]=tuningcurve(tempframes(tempframes(:,6)==0,4),...
    tempframes(tempframes(:,6)==1,4),data.samplerate);
end

%         if shuff_pass(end)==1 && plotting==1
%             figure
%             subplot(1,3,1)
%             plot(data_video_nospk(:,2),data_video_nospk(:,3),'.k');hold on
%             scatter(data_video_spk(data_video_spk(:,6)==1,2),data_video_spk(data_video_spk(:,6)==1,3),20,'r','filled');
%             axis image; box off; axis off
%
%             subplot(1,3,2)
%             imAlpha=ones(size(ratemap));
%             imAlpha(isnan(ratemap))=0;
%             imagesc(ratemap,'AlphaData',imAlpha);
%             axis xy; axis off; hold on; box off; axis image;
%             colormap(gca,viridis(255))
%             title(sprintf('true IC: %4.2f shuff IC: %4.2f',[true_IC,prctile(cell2mat(shuffled),99)]))
%
%             subplot(1,3,3)
%             h = histogram(cell2mat(shuffled),50,'FaceColor','k');hold on;
%             plot([prctile(cell2mat(shuffled),99);prctile(cell2mat(shuffled),99)],[0;max(h.BinCounts)],'k')
%             plot([true_IC;true_IC],[0;max(h.BinCounts)],'r')
%             axis square;
%             pause(.0000001)
%         end
% function IC=SHUFF(data_video_spk,secIdx,track_length)
%
% [I,row]=ismember(data_video_spk(:,1),secIdx);
% I=find(I);
% I=[1;I;length(data_video_spk)];
%
% for s=1:length(I)
%     if s==length(I)
%         datacell{s,1}=data_video_spk(I(s):I(end),:);
%         break
%     end
%     datacell{s,1}=data_video_spk(I(s):I(s+1)-1,:);
% end
%
% tempframes=data_video_spk;
% shuff_max=(length(secIdx)-20);
% bispk=[];
% IC=[];
% for ishuff=1:200
%     shift=randi([20 round(shuff_max)]);
%     tempcell=circshift(datacell,shift);
%     for open=1:length(tempcell)
%         bispk=[bispk;tempcell{open}(:,6)];
%     end
%     tempframes(:,6)=bispk;
%     bispk=[];
%
%     [ratemap,~,~,occ,~]=bindata(tempframes(tempframes(:,6)==0,:),...
%         30,tempframes(tempframes(:,6)==1,:),'no',track_length);
%
%     IC=[IC;place_cell_analysis.SpatialInformation('ratemap',...
%         ratemap,'occupancy',occ,'n_spikes',sum(data_video_spk(:,6)))];
%
% end
% end

% function IC=SHUFF(spk_ts,data_video,track_length)
% data_video(:,end)=[];
% IC=[];
% ts=linspace(0,length(data_video)/30,length(data_video));
% for i=1:200
%     shift = ((ts(end)-20)-20).*rand(1,1) + 20;
%
%     temp_spk_ts=spk_ts+shift;
%
%     idx=temp_spk_ts>data_video(end,1);
%
%     wrap_ts=(temp_spk_ts(idx)-data_video(end,1))+data_video(1,1);
%
%     temp_spk_ts(idx)=[];
%
%     temp_spk_ts=[wrap_ts;temp_spk_ts];
%
%     % FIND INDEX FOR VELOCITY FILTER (<3cm/second for 1 frame)
%     in=contiguousframes(data_video(:,5)<3,1);
%     data_video_nospk=[data_video,in];
%
%     % INTERPOLATE SPIKES TO TIMESTAMPS, POSITION, AND VEL DATA
%     X=interp1(data_video_nospk(:,1),data_video_nospk(:,2),temp_spk_ts,'linear');
%     Y=interp1(data_video_nospk(:,1),data_video_nospk(:,3),temp_spk_ts,'linear');
%
%     A=interp1(data_video_nospk(:,1),data_video_nospk(:,4),temp_spk_ts,'linear');
%     VEL=interp1(data_video_nospk(:,1),data_video_nospk(:,5),temp_spk_ts,'linear');
%     VELidx=interp1(data_video_nospk(:,1),data_video_nospk(:,6),temp_spk_ts,'linear');
%
%     % CONCAT AND SORT
%     data_video_spk=sortrows([[temp_spk_ts X Y A VEL VELidx ones(size(temp_spk_ts,1),1)];...
%         [data_video_nospk,zeros(length(data_video_nospk),1)]],1);
%     % clear TS X Y A VEL VELidx
%
%     % VELO FILTER BASED ON INDEX CREATED ABOVE
%     data_video_nospk(logical(in),:)=[];
%     data_video_nospk(:,6)=[];
%     data_video_spk(data_video_spk(:,6)==1,:)=[];
%     data_video_spk(:,6)=[];
%
%
%     [ratemap,~,~,occ,~]=bindata(data_video_nospk,...
%         30,data_video_spk(data_video_spk(:,6)==1,:),'no',track_length);
%
%     IC=[IC;place_cell_analysis.SpatialInformation('ratemap',...
%         ratemap,'occupancy',occ,'n_spikes',sum(data_video_spk(:,6)))];
%
%
% %     imAlpha=ones(size(ratemap));
% %     imAlpha(isnan(ratemap))=0;
% %     imagesc(ratemap,'AlphaData',imAlpha);
% %     axis xy; axis off; hold on; box off; axis image;
% %     colormap(gca,viridis(255))
% %     colorbar
% %     hold on
% %     pause(.05)
% end
% end