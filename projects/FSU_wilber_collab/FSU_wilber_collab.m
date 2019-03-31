function data=FSU_wilber_collab(path)
% FSU_wilber_collab main function for linear track anaysis
%
%   Input:
%           path: path to session data*
%   Output:
%           data: data structure containing tracking/spike data,
%                   measure results, and more
%
% *Example:
%   main folder: Alz39_31917
%   contents:   .t spike files
%               .nvt video file
%               .txt event file
%
% dependencies: ALLFUNCS.m, Mat2NLxVT.mexw64(or mac equivalent), ts.m, Data.m
%
% Ryan Harvey 2019


%% move to session directory
cd(path)

%% set up id & data struct
data=ALLFUNCS.datastruct(path);
clear path

%% task parameters
data.mazetype='LinearTrack';
data.mazesize=120;

%% read & parse events
data.events=ALLFUNCS.parsewilberevents;

%% read trial ts
load('trial_ts.mat','seg_times_with_AllSegments');
data.trial_ts=seg_times_with_AllSegments;
clear seg_times_with_AllSegments

%% load lfp
[data]=ALLFUNCS.handle_LFP(data);

%% extract movement info
[ts,x,y,angles] = Nlx2MatVT([data.session_path,filesep,'VT1.nvt'],[1,1,1,1,0,0],0,1);

%% check for manual coordinate correction
if exist([data.session_path,filesep,'restrictxy.mat'],'file')
    % run the following line if tracker errors from unplugs are present
    %     ALLFUNCS.manual_trackerjumps(ts,x,y,data.events(1),data.events(2),data.session_path);
    load([data.session_path,filesep,'restrictxy.mat'],'in')
    x(in==0)=NaN;
    y(in==0)=NaN;
    clear in
end

%% calculate video sample rate
data.samplerate=ALLFUNCS.samplerate(ts);

%% fix non-detects and smooth
[xtemp,ytemp]=ALLFUNCS.FixPos(x',y',ts',round(0.1667*data.samplerate));

tempangle=wrapTo360(ALLFUNCS.fixNLXangle(angles',round(0.1667*data.samplerate)))'-90;
tempangle(tempangle<0)=tempangle(tempangle<0)+360;

data.frames=[ts',xtemp,ytemp,tempangle];

clear xtemp ytemp ts tempangle x y angles

%% straighten linear track (for if you have the track at a diag)
% templinear=data.frames(data.frames(:,1)>=data.events(1) & data.frames(:,1)<=data.events(2),:);
% data.linear_track.nonlinearFrames=templinear;
% [X,Y]=ALLFUNCS.ParameterizeToLinearTrack2(templinear(:,2),templinear(:,3));
% data.frames(data.frames(:,1)>=data.events(1) & data.frames(:,1)<=data.events(2),2:3)=[X Y];
% clear X Y templinear

%% limit frames to session
data.frames=data.frames(data.frames(:,1)>=data.events(1) & data.frames(:,1)<=data.events(2),:);


%% load spikes
files=dir('*.t');
data.cellid={files.name};
data.Spikes=ALLFUNCS.LoadSpikes(data.cellid);
clear files

%% fix time: normalizes time to the first tracker frame and converts ts to seconds
data=ALLFUNCS.FixTime(data);

%% loop through each cell
for i=1:length(data.Spikes)
    
    % get matrix of spike locations
    data.spkMatrix{i}=ALLFUNCS.get_spkmatrix(data.frames,data.Spikes{i});
    
    % get overall ratemap
    [data.overall_ratemap(i,:),~,~,data.overall_occ(i,:),data.results.coherence(i)]=...
        ALLFUNCS.bindata(data.frames,data.samplerate,data.spkMatrix{i},data.mazesize);
    
    data.results.overall_InformationContent(i)=ALLFUNCS.infocontent(data.overall_ratemap(i,:),data.overall_occ(i,:));

    % get frames with embedded spike positions & spike binary 
    data.frames_w_spk{i}=sortrows([[data.spkMatrix{i},ones(size(data.spkMatrix{i},1),1)];...
        [data.frames,zeros(length(data.frames),1)]],1);
    
    % split right and left running directions
    [right,left,data.results.DirectionalityIndex(i),data.results.Displacement(i),nlaps,data.results.rateoverlap(i),...
        data.results.fieldoverlap(i),data.results.spatialcorrelation(i),startgoodrun,stopgoodrun,laps]=...
        ALLFUNCS.RightVsLeft(data.frames,data.frames_w_spk{i},data.mazesize,data.samplerate,data);
    
    if ~isfield(data.linear_track,'lapinfo')
        data.linear_track.lapinfo.startgoodrun=startgoodrun;
        data.linear_track.lapinfo.stopgoodrun=stopgoodrun;
        data.linear_track.lapinfo.laps=laps;
    end
    data.linear_track.right{i}.dataspks=right.dataspks;
    data.linear_track.right{i}.maps=right.maps;
    data.linear_track.right{i}.laps=right.laps;
    data.linear_track.left{i}.dataspks=left.dataspks;
    data.linear_track.left{i}.maps=left.maps;
    data.linear_track.left{i}.laps=left.laps;
    
    % loop through each running direction
    for iruns=1:2
        if iruns==1
            occ=right.occ; nBinsx=right.nBinsx; nBinsy=right.nBinsy;
            data.results.Coherence(iruns,i)=right.Coherence;
            SmoothRateMap=right.SmoothRateMap;
            
            [fields]=ALLFUNCS.getPlaceFields(SmoothRateMap,'minPeakRate',1,...
                'minFieldWidth',3,'percentThreshold',.2,'maxFieldWidth',length(SmoothRateMap));
            
            data.linear_track.right{1, i}.fields=fields{1, 1};
            data_video_smoothfilttemp=right.dataspks;
            data.results.lap_perm_stability(iruns,i)=right.lap_perm_stability;
            data.results.stabilityoverlaps(iruns,i)=right.stabilityoverlaps;
            data.results.meanstability(iruns,i)=right.meanstability;
        elseif iruns==2
            occ=left.occ; nBinsx=left.nBinsx; nBinsy=left.nBinsy;
            data.results.Coherence(iruns,i)=left.Coherence;
            SmoothRateMap=left.SmoothRateMap;
            
            [fields]=ALLFUNCS.getPlaceFields(SmoothRateMap,'minPeakRate',1,...
                'minFieldWidth',3,'percentThreshold',.2,'maxFieldWidth',length(SmoothRateMap));
            
            data.linear_track.left{1, i}.fields=fields{1, 1};
            data_video_smoothfilttemp=left.dataspks;
            data.results.lap_perm_stability(iruns,i)=left.lap_perm_stability;
            data.results.stabilityoverlaps(iruns,i)=left.stabilityoverlaps;
            data.results.meanstability(iruns,i)=left.meanstability;
        end
        
        data_video_spk=data_video_smoothfilttemp;
        clear ia data_video_smoothfilttemp
        
        
        %% THETA MODULATION
        [data.results.thetaindex(iruns,i),data.results.peaktheta(iruns,i),data.autocor(iruns,:,i)]=...
            ALLFUNCS.thetamodulation(data_video_spk(data_video_spk(:,end)==1,1));
        
        %% field props
        for fie=1:length(fields{1})
            PR(fie)=fields{1}{fie}.peakFR;
            start_(fie)=fields{1}{fie}.start;
            stop_(fie)=fields{1}{fie}.stop;
            temp=rescale([1:length(SmoothRateMap), fields{1, 1}{fie}.start,fields{1, 1}{fie}.stop],0,1);
            fieldbound(fie,:)=temp(end-1:end);
        end
        
        [~,I]=max(PR);
        field=zeros(1,length(SmoothRateMap));
        field(fields{1, 1}{I}.start:fields{1, 1}{I}.stop)=1;
        FieldWidth(iruns,i)=fields{1, 1}{I}.width*(data.mazesize/length(field));
        
        clear temp PR
        %% save map
        data.ratemap(iruns,:,i)=SmoothRateMap;

        %% in field / out field FR
        Map4bordertemp=SmoothRateMap;Map4bordertemp(~field)=NaN;
        data.results.infieldFR(iruns,i)=nanmean(reshape(Map4bordertemp,[],1));
 
        Map4bordertemp=SmoothRateMap;Map4bordertemp(logical(field))=NaN;
        data.results.outfieldFR(iruns,i)=nanmean(reshape(Map4bordertemp,[],1));
        clear Map4bordertemp
        
        %% spikes in field
        rescale_x=rescale(data_video_spk(:,2),0,1);
        data.results.nspikes_infield(iruns,i)=sum(data_video_spk(rescale_x>fieldbound(1) & rescale_x<fieldbound(2),end));
        
        %% LFP ANALYSIS
        % Create name of CSC file from cell and path name
        trodeID = 1;
        
        theta_phase=data.lfp.theta_phase(trodeID,...
            data.lfp.ts>=data.events(1) &...
            data.lfp.ts<=data.events(2));
        ts=data.lfp.ts(1,...
            data.lfp.ts>=data.events(1) &...
            data.lfp.ts<=data.events(2));
        theta=data.lfp.theta(trodeID,...
            data.lfp.ts>=data.events(1) &...
            data.lfp.ts<=data.events(2));
        
        % for linear track lets look at precession for every existing field
        for f=1:size(fieldbound)
            rescale_x=rescale(data_video_spk(:,2),0,1);
            temp=data_video_spk(rescale_x>fieldbound(f,1) & rescale_x<fieldbound(f,2),:);
            [ThPrecess]=ALLFUNCS.PHprecession([ts',theta_phase'],temp(temp(:,end)==1,:),temp(temp(:,end)==0,1:3),[0 1]);
            if iruns==1
                data.linear_track.right{1, i}.fields{f}.ThPrecess=ThPrecess;
            elseif iruns==2
                data.linear_track.left{1, i}.fields{f}.ThPrecess=ThPrecess;
            end
        end
        % to retain our current measures stucture, lets find (variable I)
        % the phase precession output for the field with the highest peak
        if iruns==1
            data.results.ThPrecess=data.linear_track.right{1, i}.fields{I}.ThPrecess;
        elseif iruns==2
            data.results.ThPrecess=data.linear_track.left{1, i}.fields{I}.ThPrecess;
        end
        
        [data.results.MeanThetaFreq(iruns,i),data.results.MeanOverallPow(iruns,i)]=meanfreq(theta,1000);

        clear tetrode eegfile theta_phase ts fieldbound
        
        %% calculate information content from nsma
        data.results.InformationContent(iruns,i)=ALLFUNCS.infocontent(SmoothRateMap,occ);
        
        %% Calculate the distance from place field to wall - peaks method
        [data.results.PeakRate(iruns,i),maxI]=max(SmoothRateMap);
        data.results.Field2Wall(iruns,i)=min([nBinsx-maxI maxI-0])*(data.mazesize/length(SmoothRateMap));
        
        %% OVERALL FIRING RATE
        data.results.OverallFR(iruns,i)=(sum(data_video_spk(:,end))/sum(data_video_spk(:,end)==0))*data.samplerate;
                
        %% NUMBER OF SPIKES
        data.results.nSpikes(iruns,i)=sum(data_video_spk(:,end));
        
        %% calculate sparsity from NSMA toolbox
        R=SmoothRateMap;if ~isa(R,'cell');R={R};end;RC=cat(2,R{:});[~,nCells]=size(RC);
        data.results.sparsity(iruns,i)=sum(RC,2).^2./(nCells*sum(RC.^2,2));
        
        %% calculate percent of active bins (estimate of field size - Royer et al., 2010)
        data.results.NumbActiveBins(iruns,i)=numel(find(SmoothRateMap > max(SmoothRateMap)*0.20));
        clear R rY r_max RC log_relR nBins relR occRSHP occSUM pX nCells x_max
    end
end

% save data
cd ..
cd ..
save(['ProcessedData\',data.rat,'_',data.sessionID],'-struct','data','-v7.3')

%% plot session figures

wilberFigures(data);

FolderName = data.session_path;  
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');  
  set(FigHandle,'Position',[2 42 838 924])
  print(FigHandle,'-dpng', '-r300', [FolderName,filesep,erase(FigName,{'.t',':',' '}),'.png'])
end

end