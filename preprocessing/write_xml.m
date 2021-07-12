    % write xml
    function write_xml(path,map,Fold,FNew)
    % Wrapper for buzcode bz_MakeXMLFromProbeMaps with defaults
    
    defaults.NumberOfChannels = 1;
    defaults.SampleRate = Fold;
    defaults.BitsPerSample = 16;
    defaults.VoltageRange = 20;
    defaults.Amplification = 1000;
    defaults.LfpSampleRate = FNew;
    defaults.PointsPerWaveform = 30;
    defaults.PeakPointInWaveform = 8;
    defaults.FeaturesPerWave = 4;
    [~,basename] = fileparts(path);
    bz_MakeXMLFromProbeMaps({map},path,basename,1,defaults)
    end
    
