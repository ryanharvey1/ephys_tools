function [ SmoothRateMap,nBinsx,nBinsy,occ,Coherence] = bindata(occMatrix,sampleRate,spks_VEL,linear_track,track_length)
%bindata: Creates smoothed ratemaps of linear track and open field enviorments

if isequal(linear_track, 'yes')
    nBinsx = round(track_length/3); nBinsy = 1;
    if isempty(occMatrix)
        SmoothRateMap=zeros(nBinsy,nBinsx);
        occ=zeros(nBinsy,nBinsx);
        Coherence=NaN;
        return
    end
    MinY = min(occMatrix(:,3));
    MaxY = max(occMatrix(:,3));
    MinX = min(occMatrix(:,2));
    MaxX = max(occMatrix(:,2));
    edges{1} = linspace(MinY, MaxY, nBinsy+1);
    edges{2} = linspace(MinX, MaxX, nBinsx+1);
    
    occMatrix = [occMatrix(:,3),occMatrix(:,2)];
    Omatrix = hist3(occMatrix,'Edges',edges);
    Omatrix(2,:) = [];
    Omatrix(:,end) = [];
    occ = Omatrix/sampleRate; %puts occupancy into seconds
    
    % bin spike data
    Smatrix = hist3([spks_VEL(:,3), spks_VEL(:,2)],'Edges',edges); 
    Smatrix(2,:) = [];
    Smatrix(:,end) = [];
    
    % divide binned spikes by occupancy to get rate maps and store data
    % (can use occ+eps instead of removing 0's)
    FilledRateMatrix = Smatrix./occ;
    FilledRateMatrix(isnan(FilledRateMatrix)) = 0;
    FilledRateMatrix(isinf(FilledRateMatrix)) = 0;
%     FilledRateMatrix(occ<0.150)=0;
    
    % SMOOTH
    filtWidth = [1,5]; filtSigma = 1;
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    SmoothRateMap = nanconv(FilledRateMatrix,imageFilter, 'nanout','1d');
    
    % COHERENCE
    Coherence=corr2(FilledRateMatrix,SmoothRateMap);
    
elseif isequal(linear_track,'no') % FOR ARENA
    nBinsx = round(track_length/3); nBinsy = round(track_length/3);
    MinY = min(occMatrix(:,3));
    MaxY = max(occMatrix(:,3));
    MinX = min(occMatrix(:,2));
    MaxX = max(occMatrix(:,2));
    edges{1} = linspace(MinY, MaxY, nBinsy+1);
    edges{2} = linspace(MinX, MaxX, nBinsx+1);
    
    Omatrix = hist3([occMatrix(:,3) occMatrix(:,2)],'Edges',edges);
    
    Omatrix(end,:) = [];
    Omatrix(:,end) = [];
    occ = Omatrix/sampleRate;
    %     occ(occ<0.150)=0; % Bins with less than 150ms dropped to zero
    
    % bin spike data
    Smatrix = hist3([spks_VEL(:,3), spks_VEL(:,2)],'Edges',edges);
    Smatrix(end,:) = [];
    Smatrix(:,end) = [];
    % divide binned spikes by occupancy to get rate maps and store data
    FilledRateMatrix = Smatrix./occ;
    FilledRateMatrix(isinf(FilledRateMatrix))=0;
    
    % SMOOTH
    filtWidth = [5 5]; filtSigma = 1;
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    SmoothRateMap = nanconv(FilledRateMatrix,imageFilter, 'nanout');
    
    % COHERENCE
    nonsmooth=FilledRateMatrix;
    nonsmooth(isnan(nonsmooth))=0;
    nonsmooth(isinf(nonsmooth))=0;
    
    smoothed=SmoothRateMap;
    smoothed(isnan(smoothed))=0;
    smoothed(isinf(smoothed))=0;
    
    Coherence=corr2(nonsmooth,smoothed);
end
end

