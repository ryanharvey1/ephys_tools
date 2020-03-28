function processedpath = processedpath_from_basepath(basepath)
processedpath=strsplit(basepath,filesep);
processedpath(end-2:end)=[];
processedpath = fullfile(strjoin(processedpath(1:end-3),filesep),'ProcessedData',...
    [processedpath{end-1},'_',['S',strjoin(regexp(processedpath{end},'\d*','Match'),'')]]);
end