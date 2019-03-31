cd('F:\ClarkP30_Recordings\ProcessedData');
clear
control={'ATN09','ATN08'}; %'LB01','LB03','LB05'};
transgenic={'ATN07','ATN10'}; %'LB02','LB04','LB06'};

rats=dir('*.mat');
rats={rats.name};

idx=contains(rats,["ATN09","ATN08","ATN07","ATN10"]);

rats=rats(idx);
id=[];
DR_save=[];
pr_save=[];
r_save=[];
for i=1:size(rats,2)
    disp(rats{i})
    data=load(rats{i},'frames','events','maze_size_cm','ratemap','samplerate','Spikes','spikesID');
   
    nsessions=data.events;
    if size(nsessions,2)<4
        continue
    end
    id=[id;[cellstr(repmat(rats{i},length(data.spikesID.TetrodeNum),1)),...
        data.spikesID.TetrodeNum,cellstr(num2str(data.spikesID.CellNum))]];
    
    frames=data.frames;
    
    for ses=1:4

        framestemp=frames(frames(:,1)>nsessions(1,ses) & frames(:,1)<nsessions(2,ses),:);
        for cell=1:size(data.Spikes,1)
            if ses>size(nsessions,2)
                DR_(cell,ses)=NaN;
                r_save(cell,ses)=NaN;
                pr_save(cell,ses)=NaN;

                break
            end
            [data_video_spk,~]=createframes_w_spikebinary(data,ses,cell);
            
            [r(cell,ses),~,~,pr(cell,ses),~,tuning]=tuningcurve(data_video_spk(data_video_spk(:,6)==0,4),...
                data_video_spk(data_video_spk(:,6)==1,4),data.samplerate);

            DR_(cell,ses)=distributiveRatio(data.ratemap{cell,ses},framestemp,tuning,data.maze_size_cm(ses));
        end
    end
    r_save=[r_save;r];
    pr_save=[pr_save;pr];
    DR_save=[DR_save;DR_];
    clear DR_ pr r
end


HDkeep=id(DR_save(:,1)>.2 & r_save(:,1)>.5 & pr_save(:,1)>1,:);

for i=1:size(HDkeep,1)
    data=load(HDkeep{i,1});
    p=postprocessFigures(data,{HDkeep{i,2},str2double(HDkeep(i,3))});
end


control={'ATN09','ATN08'}; %'LB01','LB03','LB05'};
transgenic={'ATN07','ATN10'}; %'LB02','LB04','LB06'};

controls=sum(contains(HDkeep(:,1),control));
Tgs=sum(contains(HDkeep(:,1),transgenic));


%__________________________LOCAL FUNCTIONS_______________________________

function DR=distributiveRatio(ratemap,frames,hdTuning,mazesize)

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
    
    PR(i)=nansum(nansum(ratemap.*theta_map))/nansum(theta_map(:));
end

for i=1:length(hdTuning)
    tempDR(i)=abs(log((1+hdTuning(i))/(1+PR(i))));
end

DR=nansum(tempDR)/length(hdTuning);

end