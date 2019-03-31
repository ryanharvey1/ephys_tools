function [ SmoothRateMap,nBinsx,nBinsy,occ,Coherence,coherenceJenni] = bindataYoder( occMatrix,sampleRate,spks_VEL,linear_track,track_length)
%bindata Summary of this function goes here
%   Detailed explanation goes here

if isequal(linear_track, 'yes')
    nBinsx = round(track_length/2); nBinsy = 1;
%     mat4occ=unique(occMatrix,'rows');
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
    FilledRateMatrix(occ<0.150)=0;
    
    % SMOOTH
    filtWidth = [1,5]; filtSigma = 1;
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    SmoothRateMap = nanconv(FilledRateMatrix,imageFilter, 'nanout','1d');
    
    % COHERENCE
    Coherence=corr2(FilledRateMatrix,SmoothRateMap);
    
    
    %
    %     % USE FMATOOLBOX TO GET RATEMAP
    %     %     map = Map([occMatrix(:,6) (occMatrix(:,2)-min(occMatrix(:,2)))/range(occMatrix(:,2))],spks_VEL(:,6),'nBins',[nBinsx nBinsy]);
    %     map = Map([occMatrix(:,5) (occMatrix(:,2)-min(occMatrix(:,2)))/range(occMatrix(:,2))],spks_VEL(:,5),'nBins',[nBinsx nBinsy],'smooth',0,'minTime',0.150);
    %     FilledRateMatrix=map.z;
    %
    %     % SMOOTH
    %     filtWidth = [1,3]; filtSigma = 1;
    %     imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    %     SmoothRateMap = nanconv(FilledRateMatrix,imageFilter, 'nanout','1d');
    %
    %     %     SmoothRateMap=map.z;
    %     occ=map.time;
    %
    %     % COHERENCE
    %     Coherence=corr2(FilledRateMatrix,SmoothRateMap);
    
elseif isequal(linear_track,'no') % FOR ARENA
    nBinsx = round(track_length/2.44); nBinsy = round(track_length/2.44);
%     mat4occ=unique(occMatrix,'rows');
    MinY = min(occMatrix(:,3));
    MaxY = max(occMatrix(:,3));
    MinX = min(occMatrix(:,2));
    MaxX = max(occMatrix(:,2));
    edges{1} = linspace(MinY, MaxY, nBinsy+1);
    edges{2} = linspace(MinX, MaxX, nBinsx+1);
    
    Omatrix = hist3([occMatrix(:,3) occMatrix(:,2)],'Edges',edges);
    
    %     Omatrix=hist3([mat4occ(:,2) mat4occ(:,3)],[nBinsy,nBinsx]);
    Omatrix(end,:) = [];
    Omatrix(:,end) = [];
    occ = Omatrix/sampleRate;
    occ(occ<0.100)=0; % Bins with less than 150ms dropped to zero

    % find pixels with zero occupancy- these will get blanked out out in the
    % rate maps
    %     Occ1NonZeroLogical = Omatrix~=0;
    %     Occ1ZeroLogical = Omatrix == 0;
    
    % bin spike data
    Smatrix = hist3([spks_VEL(:,3), spks_VEL(:,2)],'Edges',edges);
    %     Smatrix=hist3([spks_VEL(:,2), spks_VEL(:,3)],[nBinsy nBinsx]);
    Smatrix(end,:) = [];
    Smatrix(:,end) = [];
    % divide binned spikes by occupancy to get rate maps and store data
    % (can use occ+eps instead of removing 0's)
    FilledRateMatrix = Smatrix./occ;
    FilledRateMatrix(isinf(FilledRateMatrix))=0;

    
    %     % MAKE ALL VALUES OUTSIDE ARENA NAN
    %     [rowTgM,colTgM]=find(occ); %define index of non zero cells
    %     k=boundary(rowTgM,colTgM,.005); %obtain boundary of those cells
    %     boundTgMy=rowTgM(k);
    %     boundTgMx=colTgM(k);
    %     [rowTgM1,colTgM1]=find(FilledRateMatrix>=0); %obtain index of all values of matrix
    %     inTgM=inpolygon(colTgM1,rowTgM1,boundTgMx,boundTgMy); %obtain logical index of values inside boundary
    %     FilledRateMatrix(~logical(reshape(inTgM,[],length(FilledRateMatrix))))=NaN;
    %
   
    % SMOOTH
    filtWidth = [5 5]; filtSigma = 1;
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    SmoothRateMap = nanconv(FilledRateMatrix,imageFilter, 'nanout');
    

%    % SMOOTH
%     filtWidth = [5 5]; filtSigma = 1;
%     imageFilter=fspecial('gaussian',filtWidth,filtSigma);
%     for f=1:5
%         FilledRateMatrix = nanconv(FilledRateMatrix,imageFilter, 'nanout');
%     end
%     SmoothRateMap=FilledRateMatrix;
    % COHERENCE
    nonsmooth=FilledRateMatrix;
    nonsmooth(isnan(nonsmooth))=0;
    nonsmooth(isinf(nonsmooth))=0;

    smoothed=SmoothRateMap;
    smoothed(isnan(smoothed))=0;
    smoothed(isinf(smoothed))=0;

    Coherence=corr2(nonsmooth,smoothed);
    
    if max(max(nonsmooth))>100 || length(nonsmooth>=max(max(nonsmooth))*.20)==1
        test=1;
        
    end
   coherenceJenni=tiltedcoherence(nonsmooth);

    
    % Circ Rate map converts the bins outside of your circle to NaNs for cell
    %     FilledRateMatrix  = circRateMap(FilledRateMatrix,nBinsx );
end
% smooth raw binned rate matrix using a guassian
% if isequal(linear_track, 'no')
%     filtWidth = [3 3];
%     filtSigma = 1;
%     imageFilter=fspecial('gaussian',filtWidth,filtSigma);
%     SmoothRateMap = nanconv(FilledRateMatrix,imageFilter, 'nanout');
% else
%             filtWidth = [1,3];
%             filtSigma = 1;
%             imageFilter=fspecial('gaussian',filtWidth,filtSigma);
%             SmoothRateMap = nanconv(FilledRateMatrix,imageFilter, 'nanout','1d');
% end
end

