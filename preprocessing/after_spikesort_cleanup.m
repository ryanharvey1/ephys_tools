classdef after_spikesort_cleanup
    % after_spikesort_cleanup
    % handles data that was either created with spike sort 3d or MClust 4.4
    %
    % cd to session folder and run after_spikesort_cleanup.main
    %
    % Ryan Harvey
    
    methods(Static)
        
        function main
            % after_spikesort_cleanup compiles cluster data from spike sort 3d & mclust and
            % exports .mat files that contain spikes, waveforms, and quality metrics
            %
            %   Example:
            %           cd('C:\Users\ryanh\Downloads\2016-05-11_13-44-16')
            %           after_spikesort_cleanup
            %
            %   Important:
            %           make sure to include an '_clust.ntt' after your tetrode number when you save from SS3D
            %
            % Ryan Harvey
            %
            
            if exist(fullfile(pwd,'FD'),'file')
                after_spikesort_cleanup.handle_tfiles
            elseif contains(pwd,'Sorted')
                after_spikesort_cleanup.handle_ntt
            end
        end
        
        
        function handle_ntt
            fn=dir('*.ntt');
            fn=strcat(fn(1).folder,filesep,{fn.name}');
            disp(fn)
            
            for ntt=1:length(fn)
                
                % LOAD SPIKE FILE
                disp(['Loading ',fn{ntt}])
                [Timestamps,CellNumbers,Samples]=Nlx2MatSpike(fn{ntt},[1 0 1 0 1],0,1,[]);
                
                % SAVE SPIKE TIMESTAMPS
                output=[Timestamps',CellNumbers'];
                output(output(:,2)==0,:)=[];
                disp([num2str(length(unique(output(:,2)))),' Clusters'])
                disp(['Saving ',[extractBefore(fn{ntt},'_clust.ntt'),'.mat']])
                save([extractBefore(fn{ntt},'_clust.ntt'),'.mat'],'output')
                
                % CREATE INFO VARS
                cluster_n=length(unique(CellNumbers))-1;
                confidence=nan(1,cluster_n);
                final_grades=nan(1,cluster_n);
                grades=nan(cluster_n,27);
                orig_filename=erase(fn{ntt},{'Sorted\','_clust'});
                
                FD=squeeze(max(Samples))';
                clust=unique(CellNumbers);
                clust(clust<1)=[];
                ii=1;
                for i=clust
                    % ISOLATION DISTANCE
                    try
                        grades(ii,5)=IsolationDistance(FD,find(CellNumbers==i));
                    catch
                        test=1
                    end
                    % L RATIO
                    [l_output,~]=L_Ratio(FD,find(CellNumbers==i));
                    grades(ii,1)=l_output.Lratio;
                    
                    % AVERAGE WAVE FORMS
                    waves=(mean(Samples(:,:,CellNumbers==i),3)')./100;
                    means{1,ii}=[interp1(linspace(1,150,length(waves)),waves(1,:),1:150);...
                        interp1(linspace(1,150,length(waves)),waves(2,:),1:150);...
                        interp1(linspace(1,150,length(waves)),waves(3,:),1:150);...
                        interp1(linspace(1,150,length(waves)),waves(4,:),1:150)];
                    
                    % SHORT ISI
                    ISI=diff(Timestamps(CellNumbers==i)/1000) + 1e-100;
                    grades(ii,3)=sum((ISI<3))/length(ISI);
                    
                    % N SPIKES
                    grades(ii,6)=sum(CellNumbers==i);
                    ii=ii+1;
                end
                disp(['Saving ',[extractBefore(fn{ntt},'_clust.ntt'),'_info.mat']])
                save([extractBefore(fn{ntt},'_clust.ntt'),'_info.mat'],'confidence','final_grades','grades','means','orig_filename')
            end
            clear means grades
        end
        
        
        
        function handle_tfiles
            % handle data that was spike sorted in MClust 4.4
            
            % locate waveform files
            wvfn=dir('*wv.mat');
            wvfn=strcat(wvfn(1).folder,filesep,{wvfn.name}');
            for i=1:length(wvfn)
                load(wvfn{i},'mWV')
                means_all{i}=imresize(mWV,[150,size(mWV,2)])';
            end
            
            % locate cluster quality files
            qufn=dir('*CluQual.mat');
            qufn=strcat(qufn(1).folder,filesep,{qufn.name}');
            
            % CREATE INFO VARS
            cluster_n=length(unique(qufn));
            confidence_all=nan(1,cluster_n);
            final_grades_all=nan(1,cluster_n);
            grades_all=nan(cluster_n,27);
            
            
   
            orig_filename_all=strcat(extractBefore(qufn,'TT'),'TT',extractBetween(qufn,'TT','_'),'.ntt');

            
            for i=1:length(qufn)
               load(qufn{i},'CluSep') 
               grades_all(i,6)=CluSep.nSpikes;
               grades_all(i,1)=CluSep.L_Ratio.Lratio;
               grades_all(i,5)=CluSep.IsolationDistance;
            end
            
            % locate t files
            fn=dir('*.t64');
            fn=strcat(fn(1).folder,filesep,{fn.name}');
            disp(fn)
            
            for i=1:length(fn)
                path=strsplit(fn{i},filesep);
                tetrodenum(i,1)=extractBetween(path{end},'TT','_');
            end
           tetrode=unique(tetrodenum,'stable');
            
            % load spikes
            S = LoadSpikes(fn);
            
            for i=1:length(tetrode)
                nested_tetrodes{i}=[S{ismember(tetrodenum,tetrode(i))}];
            end
            
            % put Timestamps into order and create cell numbers
            for i=1:length(nested_tetrodes)
                S=nested_tetrodes{i};
                Timestamps=[];
                CellNumbers=[];
                for ii=1:length(S)
                    t=S(ii).T;
                    
                    ISI=diff(t*1000) + 1e-100;
                    ISI_store(ii,1)=sum((ISI<3))/length(ISI);
                    
                    Timestamps=[Timestamps;t];
                    CellNumbers=[CellNumbers;repmat(ii,length(t),1)];
                end
                [Timestamps,I]=sort(Timestamps);
                CellNumbers=CellNumbers(I);
                
                % convert to microseconds 
                Timestamps=Timestamps*1000000;
                
                output=[Timestamps,CellNumbers];
                
                if ~exist('Sorted','file')
                    mkdir('Sorted');
                end
                
                % SAVE SPIKE TIMESTAMPS
                disp([num2str(length(unique(output(:,2)))),' Clusters'])
                disp(['Saving TT',tetrode{i},' ',num2str(length(unique(output(:,2)))),' Clusters']) 
                
                if i==1
                    mkdir('Sorted')
                    cd Sorted
                end
                
                save(['TT',tetrode{i},'.mat'],'output')
                
                confidence=confidence_all(ismember(tetrodenum,tetrode(i)));
                final_grades=final_grades_all(ismember(tetrodenum,tetrode(i)));
                grades=grades_all(ismember(tetrodenum,tetrode(i)),:);
                means=means_all(ismember(tetrodenum,tetrode(i)));
                orig_filename=orig_filename_all{ismember(tetrodenum,tetrode(i))};
                
                grades(:,3)=ISI_store;

                save(['TT',tetrode{i},'_info.mat'],'confidence','final_grades','grades','means','orig_filename')
                
                clear ISI_store confidence final_grades grades means orig_filename
            end
            disp('finished... go post process this session :)')
        end
    end
    
end



