function mergeSessionData(sessions_to_combine)
% mergeSessionData: merges two neuralynx sessions together
% 
% Input: cell array of session ids
%        each row should contain a different session pair
%
% Example:
% sessions_to_combine={'D:\Projects\PAE_PlaceCell\data\RH11\2016-05-16_09-46-39',...
%     'D:\Projects\PAE_PlaceCell\data\RH11\2016-05-16_18-45-45'};
% mergeSessionData(sessions_to_combine)
%
% Ryan Harvey 2019

for i=1:size(sessions_to_combine,1)
    % make new folder
    mkdir([sessions_to_combine{i,1},'_combined'])
    
    % copy over cheeta log
    copyfile(fullfile(sessions_to_combine{i,1},'CheetahLogFile.txt'),...
        fullfile([sessions_to_combine{i,1},'_combined'],'CheetahLogFile.txt')); 
    
    % copy over cheetah lost ad record
    copyfile(fullfile(sessions_to_combine{i,1},'CheetahLostADRecords.txt'),...
        fullfile([sessions_to_combine{i,1},'_combined'],'CheetahLostADRecords.txt'));    

    cd(sessions_to_combine{i,1})
    
    % combine spike files
    fn=dir('*.ntt');
    fn={fn.name}';
    for ntt=1:length(fn)
        [Timestamps, ScNumbers, CellNumbers, Features, Samples, Header] =...
            Nlx2MatSpike(fullfile(sessions_to_combine{i,1},fn{ntt}), [1 1 1 1 1], 1, 1, [] );
        
        [Timestamps2, ScNumbers2, CellNumbers2, Features2, Samples2, Header2] =...
            Nlx2MatSpike(fullfile(sessions_to_combine{i,2},fn{ntt}), [1 1 1 1 1], 1, 1, [] );
        
        Mat2NlxSpike(fullfile([sessions_to_combine{i,1},'_combined'],fn{ntt}), 0, 1, [], [1 1 1 1 1],...
            [Timestamps,Timestamps2+Timestamps(end)],...
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
            Samples2, Header2] = Nlx2MatCSC(fullfile(sessions_to_combine{i,1},fn{ntt}),[1 1 1 1 1], 1, 1, [] );
        
        Mat2NlxCSC(fullfile([sessions_to_combine{i,1},'_combined'],fn{ntt}),...
            0, 1, 1, [1 1 1 1 1 1],...
            [Timestamps,Timestamps2+Timestamps(end)],...
            [ChannelNumbers,ChannelNumbers2],...
            [SampleFrequencies,SampleFrequencies2],...
            [NumberOfValidSamples,NumberOfValidSamples2],...
            [Samples,Samples2],...
            Header);
        
        disp(fullfile([sessions_to_combine{i,1},'_combined'],fn{ntt}))
    end
    
    % combine video files
    fn=dir('*.nvt');
    fn={fn.name}';
    for ntt=1:length(fn)
        [Timestamps,X,Y,Angles,Targets,Points,Header]=Nlx2MatVT(fullfile(sessions_to_combine{i,1},fn{ntt}),...
            [1 1 1 1 1 1], 1, 1, [] );
        
        [Timestamps2,X2,Y2,Angles2,Targets2,Points2,Header2]=Nlx2MatVT(fullfile(sessions_to_combine{i,2},fn{ntt}),...
            [1 1 1 1 1 1], 1, 1, [] );
        
        Mat2NlxVT(fullfile([sessions_to_combine{i,1},'_combined'],fn{ntt}),...
            0, 1, [], [1 1 1 1 1 1],...
            [Timestamps,Timestamps2+Timestamps(end)], [X,X2], [Y,Y2], [Angles,Angles2],...
            [Targets,Targets2], [Points,Points2], Header);
        
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
            [Timestamps(1),Timestamps2(2)],...
            [EventIDs(1),EventIDs2(2)],...
            TTLs,...
            Extras,...
            EventStrings,...
            Header);
        
        disp(fullfile([sessions_to_combine{i,1},'_combined'],fn{ntt}))
    end
end