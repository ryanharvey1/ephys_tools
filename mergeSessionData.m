 
function mergeSessionData
            currentFolder = pwd;
%             fn=dir('*.ntt');
            fn=dir(currentFolder);
            fn=strcat(fn(1).folder,filesep,{fn.name}');
            disp(fn)
            
            for ntt=1:length(fn)
                
                % LOAD SPIKE FILE
                disp(['Loading ',fn{ntt}])
                
                [Timestamps,Samples]=Nlx2MatSpike(fn{ntt},[1 0 0 0 1],0,1,[]);

            end
            clear means grades
        end

[Timestamps,CellNumbers,Samples]=Nlx2MatSpike(fn{ntt},[1 0 1 0 1],0,1,[]);

      
      