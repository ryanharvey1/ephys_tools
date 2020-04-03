classdef after_spikesort_cleanup
    % after_spikesort_cleanup: run after spike sorting
    % built to bring spike sort 3d, mclust 4.4, & kilosort2/phy into a
    % single format for postprocess.m to read
    %
    % Ryan Harvey
    
    methods(Static)
        
        function main
            % after_spikesort_cleanup compiles cluster data from spike sort 3d,
            % mclust 4.4, & kilosort2/phy and exports .mat files that contain 
            % spikes times, avg waveforms, and quality metrics to be loaded 
            % within load_spikes.m
            %
            %   Example:
            %           cd('C:\Users\ryanh\Downloads\2016-05-11_13-44-16')
            %           after_spikesort_cleanup.main
            %
            %   Important:
            %           make sure to include an '_clust.ntt' after your
            %           tetrode number when you save from SS3D
            %
            % Ryan Harvey
            
            
            % MClust
            if exist(fullfile(pwd,'FD'),'file')
                
                after_spikesort_cleanup.handle_tfiles
                
            % Spikesort3D    
            elseif contains(pwd,'Sorted')
                
                after_spikesort_cleanup.handle_ntt
                
            % Kilosort2 / phy    
            elseif exist(fullfile(pwd,'spike_times.npy'),'file')
                
                after_spikesort_cleanup.handle_phy
                
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
                    
                    grades(ii,5)=IsolationDistance(FD,find(CellNumbers==i));
                    
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
                save([extractBefore(fn{ntt},'_clust.ntt'),'_info.mat'],...
                    'confidence','final_grades','grades','means','orig_filename')
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
                mWV=mWV';
                for ch=1:4
                    means_all{i}(ch,:)=interp1(1:size(mWV,2),mWV(ch,:),...
                        linspace(1,size(mWV,2),150),'spline');
                end
                
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
        
        
        function handle_phy
            % handle_phy: compliles phy output into currently used .mat format
            
            % load phy output (spike times, cluster id, & waveforms)
            % calcuate cluster quality metrics
            % save all to .mat files in currently used format
            
            tic
            
            myKsDir=pwd;
            basedir = pwd;
            
            % check to see if there is a specific kilosort folder
            ks_folder = dir('Kilosort*');
            if ~isempty(ks_folder)
                myKsDir = fullfile(ks_folder.folder,ks_folder.name);
            end
            disp(myKsDir)
            
            mkdir(fullfile(pwd,'Sorted'))
            
            
            % load kilosort data processed in phy using 'spikes' function
            sp = loadKSdir(myKsDir);
            
            % read cluster_info.tsv
            opts = delimitedTextImportOptions("NumVariables", 10);
            opts.DataLines = [2, Inf];
            opts.Delimiter = "\t";
            opts.VariableNames = ["id", "Amplitude", "ContamPct", "KSLabel",...
                "amplitude", "channel", "depth", "firing_rate", "group", "n_spikes"];
            opts.VariableTypes = ["double", "double", "double", "categorical",...
                "double", "double", "double", "double", "categorical", "double"];
            opts = setvaropts(opts, 8, "TrimNonNumeric", true);
            opts = setvaropts(opts, 8, "ThousandsSeparator", ",");
            opts = setvaropts(opts, [4, 9], "EmptyFieldRule", "auto");
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule = "read";
            clusterinfo = readtable(fullfile(myKsDir,'cluster_info.tsv'), opts);
            clear opts
            
            % locate cells to keep (1=MUA, 2=Good, 3=Unsorted)
            % (Spikes from clusters labeled "noise" have already been omitted)
            spkts=[];
            clu=[];
            for i=unique(sp.clu)'
                if sp.cgs(sp.cids==i)==2
                    spkts=[spkts;sp.st(sp.clu==i)];
                    clu=[clu;sp.clu(sp.clu==i)];
                end
            end
            
            disp([num2str(sum(ismember(clusterinfo.group,'good'))),' good units'])
            disp([num2str(sum(ismember(clusterinfo.group,'mua'))),' mua units'])
            disp([num2str(sum(ismember(clusterinfo.group,'noise'))),' noise units'])
            disp('...')
            
            % set params to extract average waveforms
            datfile='filtered.dat';
            gwfparams.dataDir = [basedir,filesep];    % KiloSort/Phy output folder
            gwfparams.myKsDir = [myKsDir,filesep];
            gwfparams.fileName = datfile;         % .dat file containing the raw
            gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
            gwfparams.nCh = sp.n_channels_dat;        % Number of channels that were streamed to disk in .dat file
            gwfparams.wfWin = [-7 24]*2;              % Number of samples before and after spiketime to include in waveform
            gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out
            
            if ~isfile(datfile)
                filter_raw_dat
            end
            %             tetrodemap=reshape(repmat(1:sp.n_channels_dat/4,4,1),sp.n_channels_dat/4*4,1);
            [clusterIDs, unitQuality, contaminationRate] = sqKilosort.maskedClusterQuality(myKsDir);
            
            % Load ts from video to get first frame
            % kilosort assumes first frame is ts zero, so we need to align 
            % the spike times based on the first recorded frame, and then
            % subtract the video sample rate from the spike times to
            % account for the difference in sample rate between spikes (~32000hz)
            % and video (~30hz)
            [ts] = Nlx2MatVT(fullfile(basedir,'VT1.nvt'),[1,0,0,0,0,0],0,1);

            % split data into respective tetrodes
            for i=0:4:sp.n_channels_dat-4
                
                % spike times
                output=[spkts(ismember(clu,clusterinfo.id(ismember(clusterinfo.channel,i:i+3)))),...
                    double(clu(ismember(clu,clusterinfo.id(ismember(clusterinfo.channel,i:i+3)))))];
                
                if isempty(output)
                    continue
                end
                
                % convert seconds to microseconds, add the first video ts, 
                % and then subtract video sample rate
                output(:,1)=(output(:,1)*10^6+ts(1))-mean(diff(ts));

                disp([num2str(length(unique(output(:,2)))),' Clusters'])
                disp(['Saving ','TT',num2str(i/4+1),'.mat'])
                save(fullfile('Sorted',['TT',num2str(i/4+1),'.mat']),'output')
                
                % average waveforms
                gwfparams.spikeTimes=ceil(spkts(ismember(clu,clusterinfo.id(ismember(clusterinfo.channel,i:i+3))))*sp.sample_rate);
                gwfparams.spikeClusters = clu(ismember(clu,clusterinfo.id(ismember(clusterinfo.channel,i:i+3))));
                
                wf = getWaveForms(gwfparams);
                
                for u=1:size(wf.waveFormsMean,1)
                    ch_num=1;
                    for ch=[i+1:i+4]
                        try
                            means{u}(ch_num,:)=squeeze(wf.waveFormsMean(u,ch,:));
                        catch
                            means{u}(ch_num,:)=zeros(1,length(gwfparams.wfWin(1):gwfparams.wfWin(2)));
                        end
                        ch_num=ch_num+1;
                    end
                end
                
                orig_filename=fullfile(basedir,['TT',num2str(i/4+1),'.ntt']);
                
                confidence=NaN(1,length(means));
                
                final_grades=confidence;
                
                ui=1;
                for u=unique(output(:,2))'
                    t=((output((output(:,2)==u),1))./10^6)*1000;
                    ISI=diff(t) + 1e-100;
                    ISI_store(ui,1)=sum((ISI<3))/length(ISI);
                    ui=ui+1;
                end
                
                grades=nan(size(unique(output(:,2)),1),27);
                
                
                grades(:,1)=NaN(1,size(unique(output(:,2)),1))';
                grades(:,3)=ISI_store;
                grades(:,5)=unitQuality(ismember(clusterIDs-1,unique(output(:,2))));
                grades(:,6)=clusterinfo.n_spikes(ismember(clusterinfo.id,unique(output(:,2))));
                
                save(fullfile(basedir,'Sorted',['TT',num2str(i/4+1),'_info.mat']),'confidence','final_grades','grades','means','orig_filename')
                
                clear means grades ISI_store
                
            end
            toc
            disp('finished... go post process this session :)')
        end
    end
    
end



