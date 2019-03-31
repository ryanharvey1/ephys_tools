function DirectionalMatrix(path,event,timestamps,varargin)
%DirectionalMatrix This function creates normalized and arranged rate matrices of 1 to 2 inputs
% 
%   [ Overall_DirectionalityIndex ] = DirectionalMatrix(tfile,event,timestamps,<options>)
%   Input: 
%           tfile:      Path name for saving figures
%           event:      Number of events (sessions) currently being analyzed
%           timestamps: 1 if multible events / 0 if only one events
%           <options>:  optional list of property-value pairs (see table below)
%           =========================================================================
%           Properties    Values
%           -------------------------------------------------------------------------
%           SmoothRateMap_Right: Matrix of smoothed rate maps *
%           SmoothRateMap_Left:  Matrix of smoothed rate maps 
%           Overall_DirectionalityIndex: n x 2 matrix containing spike numbers in in each above matrix
%           group: String of group name
%           =========================================================================
%   
%   Output: 
%           Overall_DirectionalityIndex: Value of how directional your cells are
%           plots: right vs left or population matrix
%
%   Note:
%           You can include one matrix in the option field if you like, however if you want two
%           matrices, you must include Overall_DirectionalityIndex as well. 
%
%   Ryan Harvey 4/2/17

if nargin==6
    % SET VALUES
    SmoothRateMap_Right=varargin{1};
    SmoothRateMap_Left=varargin{2};
%     Overall_DirectionalityIndex=varargin{3};
    if nargin==6; group=varargin{3}; else group=[]; end
    % NORMALIZATION
    for k=1:size(SmoothRateMap_Right,1)
        for kk=1:size(SmoothRateMap_Right,2)
            SmoothRateMap_Right_Norm(k,kk) = (SmoothRateMap_Right(k,kk)-min(SmoothRateMap_Right(k,:)))/(range(SmoothRateMap_Right(k,:)));
            SmoothRateMap_Left_Norm(k,kk) = (SmoothRateMap_Left(k,kk)-min(SmoothRateMap_Left(k,:)))/(range(SmoothRateMap_Left(k,:)));
        end
    end
    % ARRANGE BINS RIGHT
    kkk=1;
    for k=1:size(SmoothRateMap_Right_Norm,2)
        index=SmoothRateMap_Right_Norm(:,k)==1;
        for kk=1:size(SmoothRateMap_Right_Norm,1)
            if index(kk,1)==1
                SmoothRateMap_Right_arrangedR(kkk,:)=SmoothRateMap_Right_Norm(kk,:);
                SmoothRateMap_Left_arrangedR(kkk,:)=SmoothRateMap_Left_Norm(kk,:);
                kkk=kkk+1;
                continue
            end
        end
    end
    
    % ARRANGE BINS LEFT
    kkk=1;
    for k=1:size(SmoothRateMap_Left_Norm,2)
        index=SmoothRateMap_Left_Norm(:,k)==1;
        for kk=1:size(SmoothRateMap_Left_Norm,1)
            if index(kk,1)==1
                SmoothRateMap_Left_arrangedL(kkk,:)=SmoothRateMap_Left_Norm(kk,:);
                SmoothRateMap_Right_arrangedL(kkk,:)=SmoothRateMap_Right_Norm(kk,:);
                kkk=kkk+1;
                continue
            end
        end
    end
    
%         SmoothRateMap_Right_arranged=imresize(SmoothRateMap_Right_arranged,[size(SmoothRateMap_Right_arranged,1)*10,size(SmoothRateMap_Right_arranged,2)*10]);
%     SmoothRateMap_Left_arranged=imresize(SmoothRateMap_Left_arranged,[size(SmoothRateMap_Left_arranged,1)*10,size(SmoothRateMap_Left_arranged,2)*10]);

%     % Overall_DirectionalityIndex
%     Overall_DirectionalityIndex=(sum(Overall_DirectionalityIndex(:,1))-...
%         sum(Overall_DirectionalityIndex(:,2)))/(sum(Overall_DirectionalityIndex(:,1))+sum(Overall_DirectionalityIndex(:,2)));
    % PLOT LEFT VS RIGHT SMOOTHED RATEMAPS
    scrsz=get(groot,'ScreenSize');
    figure2=figure('OuterPosition',[1,scrsz(4)/2,scrsz(3)/2,scrsz(4)/2]);
    % SUB 1
    subplot(2,2,1),h = pcolor(SmoothRateMap_Right_arrangedR);
    shading interp
    colormap jet
    axis off
    hold on
    colorbar
    box off
    set(h, 'EdgeColor', 'none');
    set(gca,'YDir','reverse');
    title('Smoothed Rate Map Right');
    % SUB 2
    figure2; subplot(2,2,3),h = pcolor(SmoothRateMap_Left_arrangedR);
    shading interp
    colormap jet
    axis off
    hold on
    colorbar
    box off
    set(h, 'EdgeColor', 'none');
    set(gca,'YDir','reverse');
    title('Smoothed Rate Map Left');
    % SUB 3
    figure2; subplot(2,2,2),h = pcolor(SmoothRateMap_Left_arrangedL);
    shading interp
    colormap jet
    axis off
    hold on
    colorbar
    box off
    set(h, 'EdgeColor', 'none');
    set(gca,'YDir','reverse');
    title('Smoothed Rate Map Left');
    % SUB 4
    figure2; subplot(2,2,4),h = pcolor(SmoothRateMap_Right_arrangedL);
    shading interp
    colormap jet
    axis off
    hold on
    colorbar
    box off
    set(h, 'EdgeColor', 'none');
    set(gca,'YDir','reverse');
    title('Smoothed Rate Map Right');
    
    % SAVE
    [filepath, ~] = fileparts(path);
    if timestamps==1; % SAVE WITH S# IF EVENTS EXIST
        saveas(figure2,[filepath sprintf('S%d',event) '_SmoothedRatePlot.tiff']);
        save([filepath filesep group sprintf('S%d',event) '_Data.mat']);
    else
        saveas(figure2,[filepath filesep group '_SmoothedRatePlot.tiff']);
        save([filepath filesep group '_Data.mat']);
    end
    close all
% CREATE POPULATION MATRIX
elseif nargin==4 || nargin==5
    SmoothRateMap_Right=varargin{1};
    if nargin==5; group=varargin{2}; else group=[]; end
    % NORMALIZATION
    for k=1:size(SmoothRateMap_Right,1)
        for kk=1:size(SmoothRateMap_Right,2)
            SmoothRateMap_Right_Norm(k,kk) = (SmoothRateMap_Right(k,kk)-min(SmoothRateMap_Right(k,:)))/(range(SmoothRateMap_Right(k,:)));
        end
    end
    % ARRANGE BINS
    kkk=1;
    for k=1:size(SmoothRateMap_Right_Norm,2)
        index=SmoothRateMap_Right_Norm(:,k)==1;
        for kk=1:size(SmoothRateMap_Right_Norm,1)
            if index(kk,1)==1
                SmoothRateMap_Right_arranged(kkk,:)=SmoothRateMap_Right_Norm(kk,:);
                kkk=kkk+1;
                continue
            end
        end
    end
    
%     SmoothRateMap_Right_arranged=imresize(SmoothRateMap_Right_arranged,[size(SmoothRateMap_Right_arranged,1)*10,size(SmoothRateMap_Right_arranged,2)*10]);

    % PLOT LEFT VS RIGHT SMOOTHED RATEMAPS
    scrsz=get(groot,'ScreenSize');
    figure2=figure('OuterPosition',[1,scrsz(4)/2,scrsz(3)/2,scrsz(4)/2]); h = pcolor(SmoothRateMap_Right_arranged);
    shading interp
    colormap jet
    axis off
    hold on
    colorbar
    box off
    set(h, 'EdgeColor', 'none');
    set(gca,'YDir','reverse');
    title('Smoothed Rate Matrix');
    % SAVE
    [filepath, ~] = fileparts(path);
    if timestamps==1; % SAVE WITH S# IF EVENTS EXIST
        saveas(figure2,[filepath sprintf('S%d',event) '__PopulationRateMap.tiff']);
        save([filepath filesep group sprintf('S%d',event) '_DataPoP.mat']);
    else
        saveas(figure2,[filepath filesep group '_PopulationRateMap.tiff']);
        save([filepath filesep group '_DataPoP.mat']);
    end
    close all
end
end



