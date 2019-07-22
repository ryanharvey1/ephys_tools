% HD_Project_TuningCurve_FSU_Collab
% opens .r text files containing HD cell data and creates tuning curves
%
%
% Ryan E Harvey 2018
%
% cd to data and get file names
clear;clc;close all


path='D:\Projects\HD_decoding\Ryan_RawData_Redo';

cd(path)
files=dir( '**/*.r');
folders=unique({files.folder});

for i=1:length(folders)
    area=strsplit(folders{i},filesep);
    areas{i}=strtok(area{end},' ');
end
disp(areas)

for a=1:length(folders)
    cd(folders{a})
    filenames=dir('*.r');
    filenames={filenames.name}';
    for i=1:length(filenames)
        % OPEN TEXT FILE AND CLEAN DATA
        disp(['Running:  ',filenames{i}])
        fileID=fopen(filenames{i},'r');
        dataArray=textscan(fileID,'%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]','Delimiter',...
            '\t','TextType','string','EmptyValue',NaN,'HeaderLines',2-1,'ReturnOnError',false,'EndOfLine','\r\n');
        fclose(fileID);
        frames=[dataArray{1:end-1}];
        clear dataArray fileID
        % remove last column (all NaNs)
        frames(:,end)=[];
                
        %% Add raw data to structure
        rawframes=frames;
        % remove tracker non-detects
        rawframes(rawframes(:,2)==255 | rawframes(:,3)==255 | rawframes(:,4)==255 | rawframes(:,5)==255,:)=[];
        % create expanded frames with spikes included
        framesEXP=make_framesEXP(rawframes);

        % save data
        cell = erase(filenames{i},["-",".r"," "]);
        id=areas{a};
        data.(id).(cell).frames_w_spk=framesEXP;
        data.(id).(cell).frames=rawframes;
        
        
        %% FIX POS / HEAD ANGLE / DISTANCE BETWEEN LED and add to data structure
        % remove non-detects in the first 4 xy columns
        frames(frames(:,2)==255 | frames(:,3)==255 | frames(:,4)==255 | frames(:,5)==255,2:5)=NaN;
        % add one to elimate zeros (zeros typically represent tracker non-detects)
%         frames(:,2:5)=frames(:,2:5)+1;
        % smooth position
%         [frames(:,2),frames(:,3)]=FixPos(frames(:,2),frames(:,3),frames(:,1),round(0.1667*60));
%         [frames(:,4),frames(:,5)]=FixPos(frames(:,4),frames(:,5),frames(:,1),round(0.1667*60));
        % calculate head direction 
%         frames(:,10)=wrapTo2Pi(deg2rad(atan2d(frames(:,5)-frames(:,3),frames(:,4)-frames(:,2))));
        % calculate distance between LEDs
%         frames(:,11)=sqrt((frames(:,4)-frames(:,2)).^2+(frames(:,5)-frames(:,3)).^2);
        
        % create expanded frames with spikes included
%         framesEXP_smooth=make_framesEXP(frames);

        % save data
%         data.(id).(cell).frames_w_spk_smooth=framesEXP_smooth;
%         data.(id).(cell).frames_smooth=frames;
        
        
        %% CREATE TUNING CURVE
        % 6 degree bins
        %         da=pi/30;
        %         angBins=[da/2:da:2*pi-da/2];
        %         % Occupancy
        %         histAng=hist(frames(:,10),angBins);
        %         % Number of spikes per bin
        %         spkPerAng=hist(framesEXP(framesEXP(:,6)==1,10),angBins);
        %         % Tuning
        %         hdTuning=(spkPerAng./histAng)*60;
        %         % remove nan & inf
        %         hdTuning(isnan(hdTuning) | isinf(hdTuning))=0;
        %
        %         clear da histAng spkPerAng
        
        %% COMPUTE STATS
        %         bin_centers=movmedian(angBins,2);
        %         rlength = circ_r(bin_centers',hdTuning',deg2rad(6));
        
        %% save data
        %         id=regexprep(filenames{i},'[.r]','')
%         cell = erase(filenames{i},["-",".r"," "]);
%         
%         id=areas{a};
%         data.(id).(cell).frames_w_spk=framesEXP;
%         data.(id).(cell).frames=frames;
        %         data.(id).(cell).tuningcurve=hdTuning;
        %         data.(id).(cell).rlength=rlength;
        
        
        %% PLOT
        %         tuningfig=figure; tuningfig.Color=[1 1 1];
        %         p=plot(hdTuning,'k');
        %         xlabel('Head Angle')
        %         ylabel('Firing Rate (hz)')
        %         title(['rlength: ',num2str(rlength)])
        %         set(p,'LineWidth',3)
        %         set(gca,'box','off','LineWidth',2,'XTick',linspace(0,60,7),'XTickLabel',linspace(0,60,7)*6,'FontSize',20,'FontWeight','bold')
        %
        %         area=strsplit(folders{a},filesep);
        % %         print(tuningfig,'-dpng', '-r150',['/Users/ryanharvey/Dropbox/FSU Collab/',area{6},filesep,id ' HD Cells_Tuning Curves',filesep,areas{a},'_',(cell),'.png'])
        %         close all
        %%
    end
end

function framesEXP=make_framesEXP(frames)
% remove frames with > spike binary 1 for later recombining
spksover60hz=frames(frames(:,6)>1,:);

% remove frames with > spike binary 1 for later recombining
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
end







