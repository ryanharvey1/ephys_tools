%RightVsLeft - Split up linear track data by direction & create ratemaps + directionality index
%
%   Input:
%           datavid         Linear data to be split (Must be in specfic format)*
%           dataspks        Linear data to be split that contains interpolated spike
%           Tracklength     Length of track for a 5cm bin size (will be divided by 5 to get num of x bins)
%           sampleRate      Sample rate in Hz
%   Output:
%           right                   Structure containing the following
%                                       occ:            Binned Occupancy
%                                       nBinsx:         # of bins in the x dim
%                                       nBinsy:         # of bins in the y dim
%                                       Coherence:      Cohereence measure
%                                       SmoothRateMap:  Rate map containing firing rates from only when the animal moved right
%                                       datavid:        directionally filtered linear data
%                                       dataspks:       directionally filtered linear data that contains spikes
%
%           left                    Everything above in "right", but filtered for the opposite direction
%           DirectionalityIndex     0 to 1 measure of how directionally specific the cell is
%           Displacement            Linear displacement of fields based on trajectories in # of bins
%
% *data format (the 4th column must be head angle)
%       timestamps x y angle
%       timestamps x y angle
%       timestamps x y angle
%              ... . .   ...
%
% Ryan E Harvey 2017
%
function [right,left,DirectionalityIndex,Displacement,nlaps,rateoverlap,...
    fieldoverlap,spatialcorrelation,startgoodrun,stopgoodrun,laps]=...
    RightVsLeft(data_video_nospk,data_video_spk,track_length,sampleRate,data,event)

if isfield(data.linear_track{event},'lapinfo')
    startgoodrun=data.linear_track{event}.lapinfo.startgoodrun;
    stopgoodrun=data.linear_track{event}.lapinfo.stopgoodrun;
    laps=data.linear_track{event}.lapinfo.laps;
else
    % locate laps
%     laps=FindLapsNSMAadapted(data_video_nospk(:,1),rescale(data_video_nospk(:,2),0,track_length),35);
    % locate good laps
%     [startgoodrun,stopgoodrun,laps]=NSMAFindGoodLaps(data_video_nospk(:,1),rescale(data_video_nospk(:,2),0,360),laps);
    %
    % if <15 good laps exist, dynamically loosening the edge (.1 to .3) completeness (.1 .3) threshold criteria
    % This normally should not be an issue if your animal is making nice ballistic laps
    % and your tracker is accurate.
    % note that ideally your rat should make >10 laps each session
%     if size([laps.start_ts],2)/2<15
        if track_length<100
            posbins=40;
        else
            posbins=50;
        end
        param=[[nchoosek((1:3),1),nchoosek((1:3),1)];nchoosek((1:3),2);nchoosek((3:-1:1),2)]*0.1;
        for i=1:length(param)
            laps=FindLapsNSMAadapted(data_video_nospk(:,1),rescale(data_video_nospk(:,2),0,track_length),35);
            [~,~,laps]=NSMAFindGoodLaps(data_video_nospk(:,1),rescale(data_video_nospk(:,2),0,track_length),laps,...
                param(i,1),param(i,2),0,posbins);
            nlaps(i,1)=size([laps.start_ts],2)/2;
        end
        [~,I]=max(nlaps);
        laps=FindLapsNSMAadapted(data_video_nospk(:,1),rescale(data_video_nospk(:,2),0,track_length),35);
        [startgoodrun,stopgoodrun,laps]=NSMAFindGoodLaps(data_video_nospk(:,1),rescale(data_video_nospk(:,2),0,track_length),laps,...
            param(I,1),param(I,2),0,posbins);
%     end
end

if length([laps.direction])==1
    right.occ=nan(1,track_length/3);
    right.nBinsx=track_length/3;
    right.nBinsy=1;
    right.Coherence=NaN;
    right.SmoothRateMap=nan(1,track_length/3);
    right.dataspks=nan(1,6);
    right.lap_perm_stability=NaN;
    right.stabilityoverlaps=NaN;
    right.meanstability=NaN;
    left.occ=nan(1,track_length/3);
    left.nBinsx=track_length/3;
    left.nBinsy=1;
    left.Coherence=NaN;
    left.SmoothRateMap=nan(1,track_length/3);
    left.dataspks=nan(1,6);
    left.lap_perm_stability=NaN;
    left.stabilityoverlaps=NaN;
    left.meanstability=NaN;
    DirectionalityIndex=NaN;
    Displacement=NaN;
    nlaps=1;
    rateoverlap=NaN;
    fieldoverlap=NaN;
    spatialcorrelation=NaN;
    return
end

% number of laps
nlaps=sum([laps.direction]==laps(2).direction);

% find times where rat turns around without finishing a lap to remove later
% [startgoodrun,stopgoodrun]=FindTurnAroundsNSMAadapted(data_video_nospk,laps,2);
% find laps in spike matrix

lefts=[];rights=[];i=1;
while true
    if laps(i).direction==1
        lefts=[lefts;data_video_spk(find(data_video_spk(:,1)==laps(i).start_ts):...
            find(data_video_spk(:,1)==laps(i+1).start_ts),:)];
    elseif laps(i).direction==-1
        rights=[rights;data_video_spk(find(data_video_spk(:,1)==laps(i).start_ts):...
            find(data_video_spk(:,1)==laps(i+1).start_ts),:)];
    end
    i=i+1;
    if i>=length(laps)
        break
    end
end
% remove turn around times noted above
tempdata=[rights;lefts];
[~,sortI]=sort(tempdata(:,1));
sorted=tempdata(sortI,:);
stopgoodrun(end)=sorted(end,1);
[~,Is]=intersect(sorted(:,1),startgoodrun);
[~,Ie]=intersect(sorted(:,1),stopgoodrun);
c=size(sorted,2)+1;
for i=1:length(Is)
    sorted(Is(i):Ie(i),c)=ones(length(Is(i):Ie(i)),1);
end
newInd(sortI) = 1:length(tempdata);
unsorted=sorted(newInd,:);
rights=unsorted(1:length(rights),:);
lefts=unsorted(length(rights)+1:end,:);
rights(rights(:,7)==0,:)=[];
lefts(lefts(:,7)==0,:)=[];
rights(:,end)=[];
lefts(:,end)=[];

% put laps in right and left structs
right.datavid=rights(rights(:,6)==0,:);
right.dataspks=rights;
right.spks_VEL=rights(rights(:,6)==1,:);
right.num_spikes=size(right.spks_VEL,1);

left.datavid=lefts(lefts(:,6)==0,:);
left.dataspks=lefts;
left.spks_VEL=lefts(lefts(:,6)==1,:);
left.num_spikes=size(left.spks_VEL,1);

% BIN DATA
[right.SmoothRateMap,right.nBinsx,right.nBinsy,right.occ,right.Coherence]=bindata(right.datavid,sampleRate,right.spks_VEL,'yes',track_length);

% SAME COMMENTS AS ABOVE, BUT OPPOSITE DIRECTION
[left.SmoothRateMap,left.nBinsx,left.nBinsy,left.occ,left.Coherence]=bindata(left.datavid,sampleRate,left.spks_VEL,'yes',track_length);

% CALCULATE DIRECTIONALITY INDEX (Ravassard, 2013)
% DirectionalityIndex=(right.SmoothRateMap-left.SmoothRateMap)/(right.SmoothRateMap+left.SmoothRateMap);
DirectionalityIndex=abs(sum(right.SmoothRateMap-left.SmoothRateMap)/sum(right.SmoothRateMap+left.SmoothRateMap));

% LINEAR DISPLACEMENT (Dispalcement unit is in cm)
% for ibins=0:length(right.SmoothRateMap)-1
%     SmoothRateMap_Left=circshift(left.SmoothRateMap,[0 ibins]);
%     Correlations(ibins+1,:)=corr2(right.SmoothRateMap,SmoothRateMap_Left);
% end
% [~,Displacement]=max(Correlations); Displacement=Displacement-1;

% SPATIAL CORRELATION
spatialcorrelation=corr2(right.SmoothRateMap,left.SmoothRateMap);

% RATE OVERLAP
normratemaps=rescale([right.SmoothRateMap,left.SmoothRateMap],0,1);
norm1=normratemaps(1:length(right.SmoothRateMap));
norm2=normratemaps(length(right.SmoothRateMap)+1:end);
rateoverlap=abs(max(norm1)-max(norm2));

% FIELD OVERLAP
[~,peak1]=max(right.SmoothRateMap);
[~,peak2]=max(left.SmoothRateMap);
if peak1==peak2
    fieldoverlap=1;
else
    [field1]=findfield(right.SmoothRateMap);
    [field2]=findfield(left.SmoothRateMap);
    fieldoverlap=length(intersect(find(field1),find(field2)))/sum([field1+field2]>0);
end

% DISPLACEMENT
[~,i1]=max(right.SmoothRateMap);[~,i2]=max(left.SmoothRateMap);
Displacement=abs(i1-i2)*(track_length/length(left.SmoothRateMap));

% RIGHT
% LAP STABILITY MEASURES

i=1;
mapsL=[];
mapsR=[];
corsL=[];
corsR=[];
lapsL=[];
lapsR=[];
while true
    if laps(i).direction==1
        lapsL{i,1}=left.dataspks(find(left.dataspks(:,1)==laps(i).start_ts):find(left.dataspks(:,1)==laps(i+1).start_ts),:);
        map=binlaps(left.dataspks,lapsL{i},sampleRate,lapsL{i}(lapsL{i}(:,6)==1,:),track_length);
        corsL{i,1}=corr2(left.SmoothRateMap,map);
        mapsL{i,1}=map;
    elseif laps(i).direction==-1
        lapsR{i,1}=right.dataspks(find(right.dataspks(:,1)==laps(i).start_ts):find(right.dataspks(:,1)==laps(i+1).start_ts),:);
        map=binlaps(right.dataspks,lapsR{i},sampleRate,lapsR{i}(lapsR{i}(:,6)==1,:),track_length);
        corsR{i,1}=corr2(right.SmoothRateMap,map);
        mapsR{i,1}=map;
    end
    i=i+1;
    if i>=length(laps)
        if isempty(lapsL)
            left.cors=NaN;
            left.maps=NaN;
            left.laps=NaN;
        else
            left.cors=corsL(~cellfun(@isempty,lapsL));
            left.maps=mapsL(~cellfun(@isempty,lapsL));
            left.laps=lapsL(~cellfun(@isempty,lapsL));
        end
        if isempty(lapsR)
            right.cors=NaN;
            right.maps=NaN;
            right.laps=NaN;
        else
            right.cors=corsR(~cellfun(@isempty,lapsR));
            right.maps=mapsR(~cellfun(@isempty,lapsR));
            right.laps=lapsR(~cellfun(@isempty,lapsR));
        end
        break
    end
end

% ALL PERM CORRELATIONS OF LAPS CORRELATED TO EACHOTHER
if isempty(lapsR)
    right.lap_perm_stability=NaN;
else
    corrz=[];
    for m=1:size(right.maps,1)
        for mm=1:size(right.maps,1)
            corrz=[corrz;corr2(right.maps{m},right.maps{mm})];
        end
    end
    right.lap_perm_stability=nanmean(corrz);
end
if isempty(lapsL)
    left.lap_perm_stability=NaN;
else
    corrz=[];
    for m=1:size(left.maps,1)
        for mm=1:size(left.maps,1)
            corrz=[corrz;corr2(left.maps{m},left.maps{mm})];
        end
    end
    left.lap_perm_stability=nanmean(corrz);
end

% STABILITY OVER LAPS: SEE IF LAPS ARE SIMILIARITY SIMILAR TO THE FINAL RESULTING RATE MAP OR NOT
% + CORR: BECOME STABLE
% - CORR: LOSE STABILITY
% 0 CORR: ARE STABLE THROUGHOUT LAPS
if isempty(lapsR)
    right.stabilityoverlaps=NaN;
else
    try
        right.stabilityoverlaps=corr2([1:length(right.cors)]',[right.cors{:}]');
    catch
        right.stabilityoverlaps=NaN;
    end
end
if isempty(lapsL)
    left.stabilityoverlaps=NaN;
else
    try
        left.stabilityoverlaps=corr2([1:length(left.cors)]',[left.cors{:}]');
    catch
        left.stabilityoverlaps=NaN;
    end
end
% MEAN STABILITY
if isempty(lapsR)
    right.meanstability=NaN;
else
    right.meanstability=nanmean([right.cors{:}]);
end
if isempty(lapsL)
    left.meanstability=NaN;
else
    left.meanstability=nanmean([left.cors{:}]);
end
end

function [field]=findfield(Ratemap)
[M,I]=max(Ratemap); Thres=M*0.20;
field=zeros(1,length(Ratemap)); field(I)=1;

% FORWARD
for i=1:length(Ratemap)-I
    if Ratemap(I)~=length(Ratemap)
        if Ratemap(I+i)>Thres
            field(I+i)=1;
        elseif Ratemap(I+i)<Thres
            break
        end
    end
end

% BACKWARD
for i=1:I-1
    if Ratemap(I)~=1
        if Ratemap(I-i)>Thres
            field(I-i)=1;
        elseif Ratemap(I-i)<Thres
            break
        end
    end
end
end

% function [lap,SmoothRateMap,cors]=splitlapagain(tempmat,idx,track_length,map)
% i=1;
% while true
%     if i>length(idx)
%         lap{i,1}=tempmat(idx(i-1)+1:end,:);
%         [SmoothRateMap(i,:)]=binlaps(vertcat(lap{:,:}),lap{i,1}((lap{i,1}(:,6)==0),:),...
%             30,lap{i,1}((lap{i,1}(:,6)==1),:),track_length);
%         for m=1:size(SmoothRateMap,1)
%             cors(m,1)=corr2(map,(SmoothRateMap(m,:)));
%         end
%         cors(isnan(cors))=0;
%         break
%     end
%     if i==1
%         lap{i,1}=tempmat(1:idx(i),:);
%         [SmoothRateMap(i,:)]=binlaps(vertcat(lap{:,:}),lap{i,1}((lap{i,1}(:,6)==0),:),...
%             30,lap{i,1}((lap{i,1}(:,6)==1),:),track_length);
%     else
%         lap{i,1}=tempmat(idx(i-1)+1:idx(i),:);
%         [SmoothRateMap(i,:)]=binlaps(vertcat(lap{:,:}),lap{i,1}((lap{i,1}(:,6)==0),:),...
%             30,lap{i,1}((lap{i,1}(:,6)==1),:),track_length);
%     end
%     i=i+1;
% end
% end

function [SmoothRateMap]=binlaps(occ_overall,occMatrix,sampleRate,spks_VEL,track_length)

nBinsx = round(track_length/3); nBinsy = 1;

if isempty(occMatrix)
    SmoothRateMap=zeros(nBinsy,nBinsx);
    return
end

MinY = min(occ_overall(:,3));
MaxY = max(occ_overall(:,3));
MinX = min(occ_overall(:,2));
MaxX = max(occ_overall(:,2));
edges{1} = linspace(MinY, MaxY, nBinsy+1);
edges{2} = linspace(MinX, MaxX, nBinsx+1);

occMatrix = [occMatrix(:,3),occMatrix(:,2)];
Omatrix = hist3(occMatrix,'Edges',edges);
Omatrix(2,:) = [];
Omatrix(:,end) = [];
occ = Omatrix/sampleRate;

% bin spike data
Smatrix = hist3([spks_VEL(:,3), spks_VEL(:,2)],'Edges',edges);
Smatrix(2,:) = [];
Smatrix(:,end) = [];

% divide binned spikes by occupancy to get rate maps and store data
% (can use occ+eps instead of removing 0's)
FilledRateMatrix = Smatrix./occ;
FilledRateMatrix(isnan(FilledRateMatrix)) = 0;
FilledRateMatrix(isinf(FilledRateMatrix)) = 0;

% SMOOTH
filtWidth = [1,5]; filtSigma = 1;
imageFilter=fspecial('gaussian',filtWidth,filtSigma);
SmoothRateMap = nanconv(FilledRateMatrix,imageFilter, 'nanout','1d');
end


%                                                      ___
%                                                   ,o88888
%                                                ,o8888888'
%                          ,:o:o:oooo.        ,8O88Pd8888"
%                      ,.::.::o:ooooOoOoO. ,oO8O8Pd888'"
%                    ,.:.::o:ooOoOoOO8O8OOo.8OOPd8O8O"
%                   , ..:.::o:ooOoOOOO8OOOOo.FdO8O8"
%                  , ..:.::o:ooOoOO8O888O8O,COCOO"
%                 , . ..:.::o:ooOoOOOO8OOOOCOCO"
%                  . ..:.::o:ooOoOoOO8O8OCCCC"o
%                     . ..:.::o:ooooOoCoCCC"o:o
%                     . ..:.::o:o:,cooooCo"oo:o:
%                  `   . . ..:.:cocoooo"'o:o:::'
%                  .`   . ..::ccccoc"'o:o:o:::'
%                 :.:.    ,c:cccc"':.:.:.:.:.'
%               ..:.:"'`::::c:"'..:.:.:.:.:.'
%             ...:.'.:.::::"'    . . . . .'
%            .. . ....:."' `   .  . . ''
%          . . . ...."'
%          .. . ."'
%         .
