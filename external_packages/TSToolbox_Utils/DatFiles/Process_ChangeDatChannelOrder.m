function [returnVar,msg] = Process_ChangeDatChannelOrder(fname,nbChan,ch_old,ch_new)

% USAGE:
%     Process_ChangeDatChannelOrder(fname,nbChan,ch_old,ch_new)
%     This function changes the order of a dat file
% INPUTS:
%     fname: dat file name
%     nbChan: number of channels
%     ch_old: original channel order (start at 1)
%     ch_new: new channel order (start at 1)
% 
% Adrien Peyrache 2015

[~,ix] = sort(ch_old);
ch_new = ch_new(ix);

try
    infoFile = dir(fname);
    
    chunk = 1e6;
    nbChunks = floor(infoFile.bytes/(nbChan*chunk*2));
    warning off
    if nbChunks==0
        chunk = infoFile.bytes/(nbChan*2);
    end
    m = memmapfile(fname,'Format','int16','Repeat',chunk*nbChan,'writable',true);
    d = m.Data;
    d = reshape(d,[nbChan chunk]);
    d = d(ch_new,:);
    m.Data = d(:);
    clear d m
    
    for ix=1:nbChunks-1
    %    h=waitbar(ix/(nbChunks-1));
        m = memmapfile(fname,'Format','int16','Offset',ix*chunk*nbChan*2,'Repeat',chunk*nbChan,'writable',true);
        d = m.Data;
        d = reshape(d,[nbChan chunk]);
        d = d(ch_new,:);
        m.Data = d(:);
        clear d m
    end
    %close(h)
    
    
    newchunk = infoFile.bytes/(2*nbChan)-nbChunks*chunk;
   
    if newchunk
        m = memmapfile(fname,'Format','int16','Offset',nbChunks*chunk*nbChan*2,'Repeat',newchunk*nbChan,'writable',true);
        d = m.Data;
        d = reshape(d,[nbChan newchunk]);
        d = d(ch_new,:);
        m.Data = d(:);
     clear d m
    end
    
    warning on
    returnVar = 1;
    msg = '';
    
catch
    keyboard
    returnVar = 0;
    msg = lasterr; 
end
clear m
