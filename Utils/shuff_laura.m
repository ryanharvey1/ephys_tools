function [shuff_pass,p_val,z] = shuff(groupid,varargin)
% shuff: general shuffler
%
% Input:
%       groupid: cell array with same format as below (postprocessed .mat file,
%                   tetrode id, & cell number, baseline index, cell index)
%
%     {'RH13_S20160808101427.mat','TT2.mat','2','1','4';
%     'RH13_S20160808101427.mat','TT2.mat','4','1','5';
%     'RH13_S20160808101427.mat','TT4.mat','5','1','6';
%     'RH13_S20160808101427.mat','TT5.mat','2','1','7';
%     'RH13_S20160808101427.mat','TT8.mat','3','1','8';
%     'RH13_S20160808103145.mat','TT4.mat','1','2','9';}
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
%       p_val: observed p-value
%       z: observed z-value
% ryan harvey 2019

% Output:
p = inputParser;
p.addParameter('feature',{'spatial_ic'});
p.addParameter('percentile',99);
p.addParameter('session',1);
p.addParameter('plotting',0);
p.addParameter('nshuffle',200);
p.addParameter('runningdir',[]);

p.parse(varargin{:});

feature = p.Results.feature;
percentile = p.Results.percentile;
session = p.Results.session;
nshuffle = p.Results.nshuffle;
runningdir = p.Results.runningdir;

sessions=groupid;
shuff_pass=[];

fields_to_load={'frames','events','Spikes','spikesID',...
    'maze_size_cm','sessionID','samplerate'};

if ~isempty(runningdir)
    fields_to_load={'frames','events','Spikes','spikesID',...
        'maze_size_cm','sessionID','samplerate','linear_track'};
    dirs={'right','left'};

end

for i=1:size(sessions,1)
    
    % Load session data
    data=load(sessions{i,1},fields_to_load{:});
    
    % Get cell indcies
    cells = find_cells(data,...
        sscanf(sessions{i,2},'TT%d.mat'),sessions{i,3});
    
    % Grab frames with embedded spikes for non-track
    if isempty(runningdir)
        [data_video_spk,data_video_nospk]=createframes_w_spikebinary(data,sessions{i,4},cells);
    % or track sessions
    elseif ~isempty(runningdir)
        data_video_spk=data.linear_track{1,1}.(dirs{runningdir(i)}){1,cells}.dataspks;
        data_video_nospk=data_video_spk(data_video_spk(:,6)==0,:);
    end
    
    % Evalute selected features (default spatial information
    % content) (obtain observed value)
    for f=1:length(feature)
        truth{f,1}=feval(feature{f},data_video_spk,data,sessions{i,4});
    end
    
    % Get second index
    secIdx=data_video_nospk(data.samplerate:data.samplerate:length(data_video_nospk),1);
    
    % Finds the second index within the timestamps (embedded with spikes)
    [I,~] = ismember(data_video_spk(:,1),secIdx);
    I = find(I);
    I = [1;I;length(data_video_spk)];
    
    % Breaking data up into seconds
    for s = 1:length(I)
        if s == length(I)
            datacell{s,1} = data_video_spk(I(s):I(end),:);
            break
        end
        datacell{s,1}=data_video_spk(I(s):I(s+1)-1,:);
    end
    
    % Perform shuffling
    
    % First re-order spikes by applying circular shift (random start 20 to
    % -20 seconds from end)
    tempframes=data_video_spk;
    shuff_max=(length(secIdx)-20);
    bispk=[];
    for ishuff=1:nshuffle
        shift=randi([20 round(shuff_max)]);
        tempcell=circshift(datacell,shift); % shift timestamps 
        for open=1:length(tempcell)
            bispk=[bispk;tempcell{open}(:,6)];
        end
        tempframes(:,6)=bispk;
        bispk=[];
        
        % Calculate feature(s)
        for f=1:length(feature)
            shuffled{f,ishuff}=feval(feature{f},tempframes,data,sessions{i,4});
        end
    end
    
    clear datacell
    
    % Calculate observed p-value & z-score 
    for f = 1:length(feature)
        null = cell2mat(shuffled(f,:));
        obs = truth{f,1};
        p_val(i,f) = (sum(abs(null) >= abs(obs)) + 1) / (length(null) + 1); 
        z(i,f) = (abs(obs) - abs(mean(null))) / abs(std(null));
    end
    
    shuff_pass(i,:) = cell2mat(truth)' > prctile(cell2mat(shuffled)',percentile);
    
    disp([num2str((length(shuff_pass)/length(groupid(:,1)))*100) ' percent done'])
end
end
% end

function IC=spatial_ic(tempframes,data,session)
if isfield(data,'linear_track')
    [ratemap,~,~,occ,~]=bindata(tempframes(tempframes(:,6)==0,:),...
        data.samplerate,tempframes(tempframes(:,6)==1,:),1,data.maze_size_cm(session));
else
    [ratemap,~,~,occ,~]=bindata(tempframes(tempframes(:,6)==0,:),...
        data.samplerate,tempframes(tempframes(:,6)==1,:),0,data.maze_size_cm(session));
end
IC=place_cell_analysis.SpatialInformation('ratemap',...
    ratemap,'occupancy',occ,'n_spikes',sum(tempframes(:,6)));
end
function r = mvl(tempframes,data,session)

[r,I,Ispk,peakrate,prefdirec,hdTuning]=tuningcurve(tempframes(tempframes(:,6)==0,4),...
    tempframes(tempframes(:,6)==1,4),data.samplerate);
end
function Ispk = dic(tempframes,data,session)

[r,I,Ispk,~,prefdirec,hdTuning]=tuningcurve(tempframes(tempframes(:,6)==0,4),...
    tempframes(tempframes(:,6)==1,4),data.samplerate);
end

