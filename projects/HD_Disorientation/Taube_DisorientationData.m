
%Data to play around with the Lisbon Excercise provided by Taube MIND
%Summer 2019

%% Import the data
[lr4_cell_13, ~, ~] = xlsread('d:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\HD_Disorientation\lisbon_data.xlsx','lr4 cell 13');
[lr4_cell_15, ~, ~] = xlsread('d:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\HD_Disorientation\lisbon_data.xlsx','lr4 cell 15');
[lr6_cell_5, ~, ~] = xlsread('d:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\HD_Disorientation\lisbon_data.xlsx','lr6 cell 5');
[lr4_cell_20, ~, ~] = xlsread('d:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\HD_Disorientation\lisbon_data.xlsx','lr4 cells 20-21');

cond={'Control','Dark-free','Spin 1','Spin 2','Spin 3','Spin 4'};	
subject={'lr4 cell 13','lr4 cell 15','lr6 cell 5'};
%% Organize data structure
for cell=1:3
    
    if cell==1
        temp_deg=[lr4_cell_13(1:4800,1);lr4_cell_13(1:2400,3);lr4_cell_13(1:600,5);...
            lr4_cell_13(1:600,7);lr4_cell_13(1:600,9);lr4_cell_13(1:600,11)];
        temp_spike=[lr4_cell_13(1:4800,2);lr4_cell_13(1:2400,4);lr4_cell_13(1:600,6);...
            lr4_cell_13(1:600,8);lr4_cell_13(1:600,10);lr4_cell_13(1:600,12)];
        
    elseif cell==2
        temp_deg=[lr4_cell_15(1:4800,1);lr4_cell_15(1:2400,3);lr4_cell_15(1:600,5);...
            lr4_cell_15(1:600,7);lr4_cell_15(1:600,9);lr4_cell_15(1:600,11)];
        temp_spike=[lr4_cell_15(1:4800,2);lr4_cell_15(1:2400,4);lr4_cell_15(1:600,6);...
            lr4_cell_15(1:600,8);lr4_cell_15(1:600,10);lr4_cell_15(1:600,12)];
        
    elseif cell==3
        temp_deg=[lr6_cell_5(1:4800,1);lr6_cell_5(1:2400,3);lr6_cell_5(1:600,5);...
            lr6_cell_5(1:600,7);lr6_cell_5(1:600,9);lr6_cell_5(1:600,11)];
        temp_spike=[lr6_cell_5(1:4800,2);lr6_cell_5(1:2400,4);lr6_cell_5(1:600,6);...
            lr6_cell_5(1:600,8);lr6_cell_5(1:600,10);lr6_cell_5(1:600,12)];
    end
    %% create timestamps
    frames(:,1)=linspace(0,(length(temp_spike(:,1))/10),length(temp_spike(:,1)));
    frames(:,2)=nan(size(temp_deg,1),1);
    frames(:,3)=nan(size(temp_deg,1),1);
    frames(:,4)=temp_deg;
    frames(:,5)=nan(size(temp_deg,1),1);
    frames(:,6)=temp_spike;
    
    
    events=[frames(1,1) frames(4801,1) frames(7201,1) frames(7801,1) frames(8401,1) frames(9001,1);...
        frames(4800,1) frames(7200,1) frames(7800,1) frames(8400,1) frames(9000,1) frames(end,1)];
    
    %% expand frames to create spike binary
    spksover60hz=frames(frames(:,6)>1,:);
    
    idx=find(frames(:,6)>1);
    ts=[];
    for ii=1:size(spksover60hz,1)
        if idx(ii)+1>length(frames)
            tstemp=linspace(frames(idx(ii)-1,1),frames(idx(ii),1),frames(idx(ii),6)+2);
        elseif idx(ii)-1<length(frames)
            tstemp=linspace(frames(idx(ii),1),frames(idx(ii)+1,1),frames(idx(ii),6)+2);
        else
            tstemp=linspace(frames(idx(ii)-1,1),frames(idx(ii)+1,1),frames(idx(ii),6)+2);
        end
        tstemp(1)=[];
        tstemp(end)=[];
        ts=[ts;tstemp'];
    end
    clear tstemp idx
    ts=sort([frames(frames(:,6)<=1,1);ts]);
    
    %% remove frames with > spike binary 1 for later recombining
    tempframes=frames;
    tempframes(tempframes(:,6)>1,:)=[];
    
    EXP=zeros(1,size(spksover60hz,2));
    for ii=1:size(spksover60hz,1)
        EXP=[EXP;repmat(spksover60hz(ii,:),spksover60hz(ii,6),1)];
    end
    EXP(1,:)=[];
    EXP(:,6)=ones(size(EXP,1),1);
    
    framesEXP=[tempframes;EXP];
    
    [~,I]=sort(framesEXP(:,1));
    
    framesEXP=framesEXP(I,:);
    
    framesEXP(:,1)=ts;
    

    %% Make Spike on head angle plots
    
    fig=figure;
    fig.Color=[1 1 1];
  
    for ii=1:size(events,2)
        temp_frames = framesEXP(framesEXP(:,1)>=events(1,ii) & framesEXP(:,1)<=events(2,ii),:);
        temp_rawframes = frames(frames(:,1)>=events(1,ii) & frames(:,1)<=events(2,ii),:);
        [r,~,Ispk,peakrate,prefdirec,~] = tuningcurve(temp_rawframes(:,4),temp_frames(temp_frames(:,6)==1,4),10);

        %Compute stability
        [within_Coeff,~,~] = within_HDstability(temp_frames,10);
        
        subplot(size(events,2),1,ii)
        plot(temp_frames(:,1),temp_frames(:,4),'Color',[.5 .5 .5])
        hold on;
        scatter(temp_frames(temp_frames(:,6)==1,1),temp_frames(temp_frames(:,6)==1,4),'*r')
        ylabel('head angle (deg)')
        title([subject{cell},' ',cond{ii},' stability=',num2str(within_Coeff),' ','PD=',num2str(prefdirec),' ','DirIC=',num2str(Ispk),...
            ' ','r=',num2str(r)])
        
    end
    xlabel('time (s)')

    clear temp_frames temp_rawframes
    fig=figure;
    fig.Color=[1 1 1];
    
    for ii=1:size(events,2)
        temp_frames = framesEXP(framesEXP(:,1)>=events(1,ii) & framesEXP(:,1)<=events(2,ii),:);
        temp_rawframes = frames(frames(:,1)>=events(1,ii) & frames(:,1)<=events(2,ii),:);
        [r,~,Ispk,peakrate,prefdirec,hdTuning] = tuningcurve(temp_rawframes(:,4),temp_frames(temp_frames(:,6)==1,4),10);
        subplot(size(events,2),1,ii)
        plot(6:6:360,hdTuning,'r','linewidth',2);
        ylabel('Firing Rate (Hz)')
        title([subject{cell},' ',cond{ii}])
    end
    xlabel('Head Angle (deg)')

end

clearvars -except lr4_cell_20 cond 
% observe both cells
lr4_cell_20(1,:)=[];
subject={'lr4 cell 20','lr4 cell 21'};
%% Organize data structure
for cell=1:2
    
    if cell==1
        temp_deg=[lr4_cell_20(1:4800,1);lr4_cell_20(1:2400,4);lr4_cell_20(1:600,7);...
            lr4_cell_20(1:600,10);lr4_cell_20(1:600,13);lr4_cell_20(1:600,16)];
        temp_spike=[lr4_cell_20(1:4800,2);lr4_cell_20(1:2400,5);lr4_cell_20(1:600,8);...
            lr4_cell_20(1:600,11);lr4_cell_20(1:600,14);lr4_cell_20(1:600,17)];
    elseif cell==2
        temp_deg=[lr4_cell_20(1:4800,1);lr4_cell_20(1:2400,4);lr4_cell_20(1:600,7);...
            lr4_cell_20(1:600,10);lr4_cell_20(1:600,13);lr4_cell_20(1:600,16)];
        temp_spike=[lr4_cell_20(1:4800,3);lr4_cell_20(1:2400,6);lr4_cell_20(1:600,9);...
            lr4_cell_20(1:600,12);lr4_cell_20(1:600,15);lr4_cell_20(1:600,18)];
    end
    %% create timestamps
    frames(:,1)=linspace(0,(length(temp_spike(:,1))/10),length(temp_spike(:,1)));
    frames(:,2)=nan(size(temp_deg,1),1);
    frames(:,3)=nan(size(temp_deg,1),1);
    frames(:,4)=temp_deg;
    frames(:,5)=nan(size(temp_deg,1),1);
    frames(:,6)=temp_spike;
    
    events=[frames(1,1) frames(4801,1) frames(7201,1) frames(7801,1) frames(8401,1) frames(9001,1);...
        frames(4800,1) frames(7200,1) frames(7800,1) frames(8400,1) frames(9000,1) frames(end,1)];
    %% expand frames to create spike binary
    spksover60hz=frames(frames(:,6)>1,:);
    
    idx=find(frames(:,6)>1);
    ts=[];
    for ii=1:size(spksover60hz,1)
        if idx(ii)+1>length(frames)
            tstemp=linspace(frames(idx(ii)-1,1),frames(idx(ii),1),frames(idx(ii),6)+2);
        elseif idx(ii)-1<length(frames)
            tstemp=linspace(frames(idx(ii),1),frames(idx(ii)+1,1),frames(idx(ii),6)+2);
        else
            tstemp=linspace(frames(idx(ii)-1,1),frames(idx(ii)+1,1),frames(idx(ii),6)+2);
        end
        tstemp(1)=[];
        tstemp(end)=[];
        ts=[ts;tstemp'];
    end
    clear tstemp idx
    ts=sort([frames(frames(:,6)<=1,1);ts]);
    
    %% remove frames with > spike binary 1 for later recombining
    tempframes=frames;
    tempframes(tempframes(:,6)>1,:)=[];
    
    EXP=zeros(1,size(spksover60hz,2));
    for ii=1:size(spksover60hz,1)
        EXP=[EXP;repmat(spksover60hz(ii,:),spksover60hz(ii,6),1)];
    end
    EXP(1,:)=[];
    EXP(:,6)=ones(size(EXP,1),1);
    
    framesEXP=[tempframes;EXP];
    
    [~,I]=sort(framesEXP(:,1));
    
    framesEXP=framesEXP(I,:);
    
    framesEXP(:,1)=ts;
    
    %% Make Spike on head angle plots
    
    fig=figure;
    fig.Color=[1 1 1];
  
    for ii=1:size(events,2)
        temp_frames = framesEXP(framesEXP(:,1)>=events(1,ii) & framesEXP(:,1)<=events(2,ii),:);
        temp_rawframes = frames(frames(:,1)>=events(1,ii) & frames(:,1)<=events(2,ii),:);
        [r,~,Ispk,peakrate,prefdirec,~] = tuningcurve(temp_rawframes(:,4),temp_frames(temp_frames(:,6)==1,4),10);

        %Compute stability
        [within_Coeff,~,~] = within_HDstability(temp_frames,10);
        subplot(size(events,2),1,ii)
        plot(temp_frames(:,1),temp_frames(:,4),'Color',[.5 .5 .5])
        hold on;
        scatter(temp_frames(temp_frames(:,6)==1,1),temp_frames(temp_frames(:,6)==1,4),'*r')
        ylabel('head angle (deg)')
        title([subject{cell},' ',cond{ii},' stability=',num2str(within_Coeff),' ','PD=',num2str(prefdirec),' ','DirIC=',num2str(Ispk),...
            ' ','r=',num2str(r)])
        
    end
    xlabel('time (s)')
    fig=figure;
    fig.Color=[1 1 1];
    
    for ii=1:size(events,2)
        temp_frames = framesEXP(framesEXP(:,1)>=events(1,ii) & framesEXP(:,1)<=events(2,ii),:);
        temp_rawframes = frames(frames(:,1)>=events(1,ii) & frames(:,1)<=events(2,ii),:);
        [r,~,Ispk,peakrate,prefdirec,hdTuning] = tuningcurve(temp_rawframes(:,4),temp_frames(temp_frames(:,6)==1,4),10);subaxis(size(events,2),1,ii)
        plot(6:6:360,hdTuning,'r','linewidth',2);
        title([subject{cell},' ',cond{ii}])
        ylabel('Figiring Rate (Hz)')
    end
    xlabel('Head Angle (deg)')
    
end