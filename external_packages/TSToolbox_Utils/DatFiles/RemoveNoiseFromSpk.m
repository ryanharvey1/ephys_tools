function [returnVar,msg] = RemoveNoiseFromSpk(fbasename,shankIx)

% USAGE:
%     RemoveDCfromDat(fbasename,shankIx,chanIx)
%     This function removes DC from dat files by computing the average of
%     the first 1e6 samples (or less if file is smaller)
% INPUTS:
%     fbasename: file base name
%     shankIx: vectors of shank indices
% 
% Adrien Peyrache 2011

msg = '';
returnVar = 0;

button = questdlg('Warning, this script will overwrite your SPK and RES files (back-up copies will be created), Continue?');
if strcmp(button,'Yes')
try
    readOld = '';
    shankIx = shankIx(:)';
    for shIx=shankIx
        fprintf('Shank: %i\n',shIx)
        

        fname = [fbasename '.clu.' num2str(shIx)];        
        copyF = 1;

        if exist([fname '_old'],'file')
            button = questdlg('Old files already exist, Overwrite?');
            if strcmp(button,'No')
                copyF=0;                
                button = questdlg('Load data from old files?');
                if strcmp(button,'Yes')
                    readOld = '_old';
                end
            elseif strcmp(button,'Cancel')
                m = ['\nShank ' num2str(shIx) ' not processed (canceled by user) \n'];
                fprintf(m);
                msg = [msg m];
                copyF=0;
                break;
            end
        end

        if copyF
            fprintf('     Copy files...'),tic            
            fname_old = [fname '_old'];
            eval(['!cp ' fname ' ' fname_old]);

            fname = [fbasename '.res.' num2str(shIx)];
            fname_old = [fname '_old'];
            eval(['!cp ' fname ' ' fname_old]);

            fname = [fbasename '.fet.' num2str(shIx)];
            fname_old = [fname '_old'];
            eval(['!cp ' fname ' ' fname_old]);

            fname = [fbasename '.spk.' num2str(shIx)];
            fname_old = [fname '_old'];
            eval(['!cp ' fname ' ' fname_old]);       
            fprintf('     done in %f sec. \n',toc)                 
        end
        
        fname = [fbasename '.clu.' num2str(shIx)];      
        clu = load(fname);
        nClu = clu(1);
        clu = clu(2:end);
        cluIx = clu>0;      
        
        %Write new bad file
        fprintf('     Bad Ix file...'),tic
        badIx = clu==0;
        fname = [fbasename '.bad.' num2str(shIx)];                
        fid = fopen([fname],'w');       
        fprintf(fid,'%i\n',find(badIx));
        fclose(fid);

        fprintf('  done in %f sec. \n',toc)                 
       
        %Write new clu file (unnecassry, but useful to test the prog with
        %kluster)
        fprintf('     New Clu file...'),tic
        fname = [fbasename '.clu.' num2str(shIx)];                
        fid = fopen([fname '.tmp'],'w');  
        fprintf(fid,'%i\n',nClu-1);        
        fprintf(fid,'%i\n',clu(cluIx));
        fclose(fid);
        clear clu
        fprintf('  done in %f sec. \n',toc)   

        %Write new res file
        fprintf('     New Res file...'),tic
        fname = [fbasename '.res.' num2str(shIx)];
        res = load(fname);
        fid = fopen(fname,'w');  
        fprintf(fid,'%i\n',res(cluIx));
        fclose(fid);
        fprintf('  done in %f sec. \n',toc)   
        
        %Write new fet
        fprintf('     New Fet file...'),tic
        fname = [fbasename '.fet.' num2str(shIx)];
        fet = dlmread(fname,' ',1,0);
        fet=fet(cluIx,[1:end-2 end]);
        ch = '';
        nbDim = size(fet,2);
        for ii=1:size(fet,2)
            ch = [ch ' %i'];
        end
        fet = fet';
        fet = fet(:);
        fid = fopen([fname],'w');
        fprintf(fid,'%i\n',nbDim);
        fprintf(fid,[ch '\n'],fet);
        fclose(fid);
        clear fet
        fprintf('  done in %f sec. \n',toc)        
        
        
        %Write Spk fet
        fprintf('     New Spk file...'),tic
        fname = [fbasename '.spk.' num2str(shIx)];
        
        maxSamplesPerChunk = 2^20;
        nChannels=8;
        frequency=20000;
        sizeInBytes = 2;
        channels = 1:nChannels;
        skip = 0;
        f = fopen(fname,'r');
        if f == -1,
            error(['Cannot read ' filename ' (insufficient access rights?).']);
        end

        fileStart = ftell(f);
        status = fseek(f,0,'eof');
        if status ~= 0,
            fclose(f);
            error('Error reading the data file (possible reasons include trying to read past the end of the file).');
        end
        fileStop = ftell(f);
        nSamplesPerChannel = (fileStop-fileStart)/nChannels/sizeInBytes;
        nSamples = nChannels*nSamplesPerChannel;
        duration = nSamplesPerChannel/frequency;
        frewind(f);
        status = fseek(f,0,'bof');
        if status ~= 0,
            fclose(f);
            error('Could not start reading (possible reasons include trying to read past the end of the file).');
        end
        nSamplesPerChunk = floor(maxSamplesPerChunk/nChannels)*nChannels;
        durationPerChunk = nSamplesPerChunk/frequency/nChannels;
        nChunks = floor(duration/durationPerChunk);
        fprintf('nb chunks: %d',nChunks);
        fo = fopen([fname '.tmp'],'w');
        % Read all chunks
        i = 1;
        for j = 1:nChunks,
            if ~mod(j/nChunks,0.05)
                fprintf('.')
            end
            
            d = LoadBinaryChunk(f,'frequency',frequency,'nChannels',nChannels,'channels',channels,'duration',durationPerChunk,'skip',skip);
           
            d = reshape(d',nChannels,32,[]);
            [l,m,n] = size(d);
            ix = badIx(i:i+n-1);
            if n == 0, break; end
            d(:,:,ix) = [];
            fwrite(fo,d(:),'int16');
            i = i+n;
        end
        % If the data size is not a multiple of the chunk size, read the remainder
        remainder = duration - nChunks*durationPerChunk;
        if remainder ~= 0,
            d = fread(f,[nChannels inf],'int16');
            d = reshape(d,nChannels,32,[]);
            [l,m,n] = size(d);
            ix = badIx(i:i+n-1);
            if n == 0, break; end
            d(:,:,ix) = [];
            fwrite(fo,d(:),'int16');
        end
        fclose(fo);
        fclose(f);
        
        eval(['!mv ' fname '.tmp ' fname ';']);
        fprintf('     done in %f sec. \n',toc);
        
    end
    
catch
%     returnVar = 0;
%     msg = lasterr; 
keyboard
end
end
