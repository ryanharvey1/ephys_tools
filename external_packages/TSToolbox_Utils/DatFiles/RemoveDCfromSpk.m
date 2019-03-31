function [returnVar,msg] = RemoveDCfromSpk(fbasename,chan)

% USAGE:
%     RemoveDCfromDat(fbasename,shankIx,chanIx)
%     This function removes DC from dat files by computing the average of
%     the first 1e6 samples (or less if file is smaller)
% INPUTS:
%     fname: dat file name
%     nbChan: total number of channels in dat file
%     chanIx: vectors of channel indices
% 
% Adrien Peyrache 2011

try

    %data = memmapfile(fname,'format','int16','Writable',true);
    spk = LoadSpikeWaveforms([fbasename '.spk.' num2str(chan)],8,32);
    for ii=1:8
        m = mean(squeeze(spk(ii,:,:)));
        warning off
        spk(ii,:,:) = squeeze(spk(ii,:,:))-int16(ones(32,1)*m);
        warning on;
    end
    spk = spk(:);
    fid = fopen([fbasename '.spk.' num2str(chan)],'r');
    fwrite(fid,spk,'int16');
    fclose(fid);
   
    returnVar = 1;
    msg = '';
    
catch
    returnVar = 0;
    msg = lasterr; 
end
