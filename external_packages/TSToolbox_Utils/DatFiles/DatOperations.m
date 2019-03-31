function DatAnalysis(rootFolder,fction,fbasename,datIx,varargin)

% USAGE:
%     DatOperations(fbasename,shankIx,chanIx)
% INPUTS:
%     rootFolder: root folder of the data
%     fction: name of function to be called
%     fbasename: file base name
%     shankIx: vector of shank indices to be treated
%     chanIx (optionnal): vectors of channel indices (will be the same for all shanks)
% 
% The program assumes that the data are stored in the folder "rootFolder" 
%
% Adrien Peyrache 2011

datIx = datIx(:)';

disp(['Launch ' fction ' for ' fbasename])
for ix=datIx
    disp(['Dat file ' num2str(ix)])
    fname = [rootFolder filesep fbasename '-0' num2str(ix) '.dat'];
    
    tic
    eval(['[returnVar,msg] = ' fction '(''' fname ''',33,[1:32]);']);
    timeSpt = toc;
    
    if returnVar
        disp(['...done in ' num2str(timeSpt) ' sec.'])
    else
        disp(['...error:' msg])
    end
end
