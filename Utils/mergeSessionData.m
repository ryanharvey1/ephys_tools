function mergeSessionData(sessions_to_combine,varargin)
% mergeSessionData: merges two neuralynx sessions together
% 
% Input: cell array of session ids
%        each row should contain a different session pair
%
% optional: 
%       'view_xy': 1 to view xy coordinates from both sessions before merging. (default 0)
%  
% Example:
% sessions_to_combine={'D:\Projects\PAE_PlaceCell\data\RH11\2016-05-16_09-46-39',...
%     'D:\Projects\PAE_PlaceCell\data\RH11\2016-05-16_18-45-45'};
% mergeSessionData(sessions_to_combine)
%
% Ryan Harvey 2019
warning off

p = inputParser;
p.addParameter('view_xy',0);
p.parse(varargin{:});

view_xy = p.Results.view_xy;

for i=1:size(sessions_to_combine,1)
    
    disp([sessions_to_combine{i,1},'  &&  ',sessions_to_combine{i,2}])
    
    % compare xy from each session pair 
    if view_xy
        cd(sessions_to_combine{i,1})
        fn=dir('*.nvt');
        fn={fn.name}';
                
        [Timestamps,X,Y,Angles,Targets,Points,Header]=Nlx2MatVT(fullfile(sessions_to_combine{i,1},fn{:}),...
            [1 1 1 1 1 1], 1, 1, [] );
        
        [Timestamps2,X2,Y2,Angles2,Targets2,Points2,Header2]=Nlx2MatVT(fullfile(sessions_to_combine{i,2},fn{:}),...
            [1 1 1 1 1 1], 1, 1, [] );
        
        fig=figure;
        subplot(1,2,1)
        plot(X,Y,'.k')
        xlim([min([X,X2]) max([X,X2])])
        ylim([min([Y,Y2]) max([Y,Y2])])
        title(sessions_to_combine{i,1})
        grid on
        
        subplot(1,2,2)
        plot(X2,Y2,'.k')
        xlim([min([X,X2]) max([X,X2])])
        ylim([min([Y,Y2]) max([Y,Y2])])
        title(sessions_to_combine{i,2})
        grid on
        
        str = input('Continue to merge? y or n  ','s');
        close all
        if strcmp(str,'n')
            continue
        end
    end
    
    % make new folder
    mkdir([sessions_to_combine{i,1},'_combined'])
    
    % copy over cheeta log
    copyfile(fullfile(sessions_to_combine{i,1},'CheetahLogFile.txt'),...
        fullfile([sessions_to_combine{i,1},'_combined'],'CheetahLogFile.txt')); 
    
    % copy over cheetah lost ad record
    copyfile(fullfile(sessions_to_combine{i,1},'CheetahLostADRecords.txt'),...
        fullfile([sessions_to_combine{i,1},'_combined'],'CheetahLostADRecords.txt'));    

    cd(sessions_to_combine{i,1})
    
    % combine video files
    fn=dir('*.nvt');
    fn={fn.name}';
    for ntt=1:length(fn)
        [Timestamps,X,Y,Angles,Targets,Points,Header]=Nlx2MatVT(fullfile(sessions_to_combine{i,1},fn{ntt}),...
            [1 1 1 1 1 1], 1, 1, [] );
        
        [Timestamps2,X2,Y2,Angles2,Targets2,Points2,Header2]=Nlx2MatVT(fullfile(sessions_to_combine{i,2},fn{ntt}),...
            [1 1 1 1 1 1], 1, 1, [] );
        
        % time gap of time between video frames
        ts_gap=mean(diff(Timestamps));
        
        % range of time between sessions
        ts_range=[Timestamps(end),Timestamps2(1)]; 
        
        % alter second session time stamps to fit after first session
        tempts=((Timestamps2-ts_range(2))+ts_range(1))+ts_gap;
        
        % locate final time stamp
        lastts=tempts(end);
        
        Mat2NlxVT(fullfile([sessions_to_combine{i,1},'_combined'],fn{ntt}),...
            0, 1, [], [1 1 1 1 1 1],...
            [Timestamps,tempts], [X,X2], [Y,Y2], [Angles,Angles2],...
            [Targets,Targets2], [Points,Points2], Header);
        
        disp(fullfile([sessions_to_combine{i,1},'_combined'],fn{ntt}))
    end
    
    % combine spike files
    fn=dir('*.ntt');
    fn={fn.name}';
    for ntt=1:length(fn)
        [Timestamps, ScNumbers, CellNumbers, Features, Samples, Header] =...
            Nlx2MatSpike(fullfile(sessions_to_combine{i,1},fn{ntt}), [1 1 1 1 1], 1, 1, [] );
        
        [Timestamps2, ScNumbers2, CellNumbers2, Features2, Samples2, Header2] =...
            Nlx2MatSpike(fullfile(sessions_to_combine{i,2},fn{ntt}), [1 1 1 1 1], 1, 1, [] );
        
        tempts=(Timestamps2-ts_range(2));
        to_remove=tempts<0;
        tempts(to_remove)=[];
        tempts=tempts+ts_range(1)+ts_gap;
        
        ScNumbers2(to_remove)=[];
        CellNumbers2(to_remove)=[];
        Features2(:,to_remove)=[];
        Samples2(:,:,to_remove)=[];
        
        Mat2NlxSpike(fullfile([sessions_to_combine{i,1},'_combined'],fn{ntt}), 0, 1, [], [1 1 1 1 1],...
            [Timestamps,tempts],...
            [ScNumbers,ScNumbers2],...
            [CellNumbers,CellNumbers2],...
            [Features,Features2],...
            cat(3,Samples,Samples2),...
            Header);
        
        disp(fullfile([sessions_to_combine{i,1},'_combined'],fn{ntt}))
    end
    
    % combine csc files
    fn=dir('*.ncs');
    fn={fn.name}';
    for ntt=1:length(fn)
        [Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples,...
            Samples, Header] = Nlx2MatCSC(fullfile(sessions_to_combine{i,1},fn{ntt}),[1 1 1 1 1], 1, 1, [] );
        
        [Timestamps2, ChannelNumbers2, SampleFrequencies2, NumberOfValidSamples2,...
            Samples2, Header2] = Nlx2MatCSC(fullfile(sessions_to_combine{i,2},fn{ntt}),[1 1 1 1 1], 1, 1, [] );
        
        tempts=(Timestamps2-ts_range(2));
        to_remove=tempts<0;
        tempts(to_remove)=[];
        tempts=tempts+ts_range(1)+ts_gap;
        
        ChannelNumbers2(to_remove)=[];
        SampleFrequencies2(to_remove)=[];
        NumberOfValidSamples2(to_remove)=[];
        Samples2(:,to_remove)=[];
        
        Mat2NlxCSC(fullfile([sessions_to_combine{i,1},'_combined'],fn{ntt}),...
            0, 1, 1, [1 1 1 1 1 1],...
            [Timestamps,tempts],...
            [ChannelNumbers,ChannelNumbers2],...
            [SampleFrequencies,SampleFrequencies2],...
            [NumberOfValidSamples,NumberOfValidSamples2],...
            [Samples,Samples2],...
            Header);
        
        disp(fullfile([sessions_to_combine{i,1},'_combined'],fn{ntt}))
    end
    

    % combine event files
    fn=dir('*.nev');
    fn={fn.name}';
    for ntt=1:length(fn)
        [Timestamps,EventIDs,TTLs,Extras,EventStrings,Header]=...
            Nlx2MatEV(fullfile(sessions_to_combine{i,1},fn{ntt}), [1 1 1 1 1], 1, 1, [] );
        
        [Timestamps2,EventIDs2,TTLs2,Extras2,EventStrings2,Header2]=...
            Nlx2MatEV(fullfile(sessions_to_combine{i,1},fn{ntt}), [1 1 1 1 1], 1, 1, [] );
        
        Mat2NlxEV(fullfile([sessions_to_combine{i,1},'_combined'],fn{ntt}),...
            0, 1, [], [1 1 1 1 1],...
            [Timestamps(1),lastts],...
            [EventIDs(1),EventIDs2(2)],...
            TTLs,...
            Extras,...
            [EventStrings(1);EventStrings(2)],...
            Header);
        
        disp(fullfile([sessions_to_combine{i,1},'_combined'],fn{ntt}))
    end
end
warning on
end